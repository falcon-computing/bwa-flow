#include <algorithm>
#include <boost/asio.hpp>
#include <boost/function_types/result_type.hpp>
#include <boost/make_shared.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread.hpp>
#include <fstream>
#include <list>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "bwa/utils.h"

#ifdef USE_HTSLIB
#include "htslib/ksort.h"
#endif

#include "bwa_wrapper.h"
#include "config.h"
#include "Pipeline.h"  
#include "util.h"

// Comparator function for bam1_t records
#ifdef USE_HTSLIB
bool bam1_lt(const bam1_t* a, const bam1_t* b) {
  return ((uint64_t)a->core.tid<<32|(a->core.pos+1)<<1|bam_is_rev(a)) 
       < ((uint64_t)b->core.tid<<32|(b->core.pos+1)<<1|bam_is_rev(b));
}
#endif

#ifdef SCALE_OUT
#include "mpi.h"
// Encode a scalar value to serialized data
template <typename T>
static inline void putT(std::stringstream &ss, T value) {
  ss.write(reinterpret_cast<char*>(&value), sizeof(value));
}

// Decode a scalar value from serialized data
template <typename T>
static inline void getT(std::stringstream &ss, T &value) {
  ss.read(reinterpret_cast<char*>(&value), sizeof(value));
}

// Store a string with its length to serialized data
static inline void putStr(std::stringstream &ss, const char* str) {
  if (str) {
    size_t length = strlen(str);
    putT(ss, length);
    ss.write(str, length);
  }
  else {
    size_t length = 0;
    putT(ss, length);
  }
}

// Retrieve a string from serialized data
static inline void getStr(std::stringstream &ss, char* &dst) {
  size_t length = 0;
  getT(ss, length);

  if (length > 0) {
    dst = (char*)malloc(length+1);
    ss.read(dst, length);
    dst[length] = '\0';
  }
}

// Check MPI::Request status with mutex lock
static inline bool queryStatus(MPI::Request &req) {
  boost::mutex::scoped_lock lock(mpi_mutex);
  return req.Test();
}

static inline void bwaMPISend(const void* buf,
    int count, const MPI::Datatype& datatype,
    int dest, int tag
) {
  MPI::Request req;
  {
    boost::mutex::scoped_lock lock(mpi_mutex);
    req = MPI::COMM_WORLD.Isend(buf, count,
        datatype, dest, tag);
  }
  while (!queryStatus(req)) {
    boost::this_thread::sleep_for(boost::chrono::milliseconds(1));
  }
}

static inline void bwaMPIRecv(void* buf,
    int count, const MPI::Datatype& datatype,
    int source, int tag
) {
  MPI::Request req;
  {
    boost::mutex::scoped_lock lock(mpi_mutex);
    req = MPI::COMM_WORLD.Irecv(buf, count,
        datatype, source, tag);
  }
  while (!queryStatus(req)) {
    boost::this_thread::sleep_for(boost::chrono::milliseconds(1));
  }
}

void SeqsDispatch::compute() {

  int rank   = mpi_rank;
  int nprocs = mpi_nprocs;

  uint64_t num_seqs_produced = 0;

  bool finished = false;
  while (!finished) { 
    SeqsRecord record;
    bool ready = this->getInput(record);

    while (!this->isFinal() && !ready) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(100));
      ready = this->getInput(record);
    }
    if (ready) { // record is a valid new input
      uint64_t start_ts = getUs();

      // Serialize output record
      std::string ser_data = serialize(&record);

      int length = ser_data.length();

      if (length <= 0) {
        throw std::runtime_error("Possible overflow of msg length");
      }

      for (int i = 0; i < record.batch_num; i++) {
        free(record.seqs[i].name);
        free(record.seqs[i].comment);
        free(record.seqs[i].seq);
        free(record.seqs[i].qual);
      }
      free(record.seqs);

      DLOG_IF(INFO, VLOG_IS_ON(2)) << "Serializing seq batch of "
        << length / 1024  << "kb"
        << " in " << getUs() - start_ts << " us";

      start_ts = getUs();

      // First query a process to send data to
      int proc_id = -1;
      bwaMPIRecv(&proc_id, 1, MPI::INT, MPI::ANY_SOURCE, SEQ_DP_QUERY);

      bwaMPISend(&length, 1, MPI::INT, proc_id, SEQ_DP_LENGTH);

      bwaMPISend(ser_data.c_str(), length, MPI::CHAR, proc_id, SEQ_DP_DATA);

      DLOG_IF(INFO, VLOG_IS_ON(1)) << "Sending seqs batch " << record.start_idx
        << " to proc_" << proc_id
        << " takes " << getUs() - start_ts << " us";
    }
    else {
      // this means isFinal() is true and input queue is empty
      DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finish reading seqs, start send finish signals";

      for (int p = 0; p < nprocs; p++) {
        int proc_id = 0;
        int length = 0;
        bwaMPIRecv(&proc_id, 1, MPI::INT, MPI::ANY_SOURCE, SEQ_DP_QUERY);

        // Send finish signal to child-process
        bwaMPISend(&length, 1, MPI::INT, p, SEQ_DP_LENGTH);
      }
      finished = true;
    }
  }
}

std::string SeqsDispatch::serialize(SeqsRecord* data) {

  uint64_t start_idx = data->start_idx;
  int      batch_num = data->batch_num;

  std::stringstream ss;

  putT(ss, start_idx);
  putT(ss, batch_num);

  for (int i = 0; i < batch_num; i++) {
    bseq1_t* seq = &data->seqs[i];

    putStr(ss, seq->name);
    putStr(ss, seq->comment);
    putStr(ss, seq->seq);
    putStr(ss, seq->qual);
  }

  return ss.str();
}

void SeqsReceive::compute() {

  int rank   = mpi_rank;
  int nprocs = mpi_nprocs;

  bool finished = false;
  while (!finished) {
    uint64_t start_ts = getUs();

    // Request new data from master
    bwaMPISend(&rank, 1, MPI::INT, MASTER_RANK, SEQ_DP_QUERY);

    int length = 0;
    bwaMPIRecv(&length, 1,
        MPI::INT, MASTER_RANK, SEQ_DP_LENGTH);

    if (length > 0) {
      char* ser_data = (char*) malloc(length);

      bwaMPIRecv(ser_data, length,
          MPI::CHAR, MASTER_RANK, SEQ_DP_DATA);

      SeqsRecord output = deserialize(ser_data, length);
      free(ser_data);

      DLOG_IF(INFO, VLOG_IS_ON(1)) << "Receive one read batch in "
        << getUs() - start_ts << " us";

      pushOutput(output);
    }
    else {
      // Means master has no more batch
      finished = true;
    }
  }
}

SeqsRecord SeqsReceive::deserialize(const char* data, size_t length) {

  uint64_t start_idx = 0;
  int      batch_num = 0;

  std::stringstream ss;
  ss.write(data, length);

  getT(ss, start_idx);
  getT(ss, batch_num);

  bseq1_t* seqs = (bseq1_t*)malloc(batch_num*sizeof(bseq1_t));

  for (int i = 0; i < batch_num; i++) {
    bseq1_t* seq = &seqs[i];
    memset(seq, 0, sizeof(bseq1_t));

    getStr(ss, seq->name);
    getStr(ss, seq->comment);
    getStr(ss, seq->seq);
    getStr(ss, seq->qual);

    seq->l_seq = strlen(seq->seq);
  }

  SeqsRecord output;
  output.start_idx = start_idx;
  output.batch_num = batch_num;
  output.seqs = seqs;

  return output;
}

void SamsSend::compute() {

#ifndef USE_HTSLIB    
  int rank   = mpi_rank;
  int nprocs = mpi_nprocs;

  bool finished = false;
  while (!finished) {
    SeqsRecord input;
    bool ready = this->getInput(input);

    while (!this->isFinal() && !ready) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(100));
      ready = this->getInput(input);
    }

    if (!ready) { 
      // This means isFinal() is true and input queue is empty
      int length = 0;

      // Send proc_id to master to let master receive following msg
      bwaMPISend(&rank, 1, MPI::INT, MASTER_RANK, SAM_RV_QUERY);

      // Send a zero-length to master indicating current process is finished
      bwaMPISend(&length, 1, MPI::INT, MASTER_RANK, SAM_RV_LENGTH);

      finished = true;
    }
    else {
      uint64_t start_ts = getUs();

      // Serialize data and send to master
      std::string ser_data = serialize(&input);

      int length = ser_data.length();

      if (length <= 0) {
        throw std::runtime_error("Possible overflow of msg length");
      }
      DLOG_IF(INFO, VLOG_IS_ON(2)) << "Serializing sam batch of "
        << length / 1024  << "kb"
        << " in " << getUs() - start_ts << " us";

      start_ts = getUs();

      // Send proc_id to master to let master receive following msg
      bwaMPISend(&rank, 1, MPI::INT, MASTER_RANK, SAM_RV_QUERY);

      bwaMPISend(&length, 1, MPI::INT, MASTER_RANK, SAM_RV_LENGTH);

      bwaMPISend(ser_data.c_str(), length,
          MPI::CHAR, MASTER_RANK, SAM_RV_DATA);

      DLOG_IF(INFO, VLOG_IS_ON(1)) << "Sending sam batch " << input.start_idx
        << " to master takes " << getUs() - start_ts << " us";

      //freeSeqs(input.seqs, input.batch_num);
      for (int i = 0; i < input.batch_num; i++) {
        free(input.seqs[i].sam);
      }
      free(input.seqs);
    }
  }
#else
  LOG(ERROR) << "SamSend() is not supported in sorted bwa version";
#endif
}

std::string SamsSend::serialize(SeqsRecord* data) {

#ifndef USE_HTSLIB    
  uint64_t start_idx = data->start_idx;
  int      batch_num = data->batch_num;
  
  std::stringstream ss;

  putT(ss, start_idx);
  putT(ss, batch_num);

  for (int i = 0; i < batch_num; i++) {
    bseq1_t* seq = &data->seqs[i];

    putStr(ss, seq->sam);
  }

  return ss.str();
#endif
}

void SamsReceive::compute() {
  
#ifndef USE_HTSLIB    
  int rank   = mpi_rank;
  int nprocs = mpi_nprocs;

  // Recording finish status of each child process
  std::unordered_map<int, bool> unfinished_proc;

  for (int p = 0; p < nprocs; p++) {
    unfinished_proc[p] = true;
  }

  while (!unfinished_proc.empty()) {

    int proc_id = -1;
    int length = 0;

    // non-blocking query for tasks
    bwaMPIRecv(&proc_id, 1, MPI::INT, MPI::ANY_SOURCE, SAM_RV_QUERY);

    bwaMPIRecv(&length, 1, MPI::INT, proc_id, SAM_RV_LENGTH);

    if (length == 0) {
      // Process proc_id is already finished, remove from table
      unfinished_proc.erase(proc_id);
    }
    else {
      // Allocate buffer for serialized obj
      char* ser_data = (char*) malloc(length);

      bwaMPIRecv(ser_data, length, MPI::CHAR, proc_id, SAM_RV_DATA);

      SeqsRecord output = deserialize(ser_data, length);
      free(ser_data);

      pushOutput(output);
    }
  }
#else
  LOG(ERROR) << "SamSend() is not supported in sorted bwa version";
#endif
}

SeqsRecord SamsReceive::deserialize(const char* data, size_t length) {

#ifndef USE_HTSLIB    
  uint64_t start_idx = 0;
  int      batch_num = 0;

  std::stringstream ss;
  ss.write(data, length);
  
  // Parse integers from serialized data
  getT(ss, start_idx);
  getT(ss, batch_num);

  bseq1_t* seqs = (bseq1_t*)malloc(batch_num*sizeof(bseq1_t));

  for (int i = 0; i < batch_num; i++) {
    bseq1_t* seq = &seqs[i];
    memset(seq, 0, sizeof(bseq1_t));
    
    // Parse one string from serialized data
    getStr(ss, seq->sam);
  }

  SeqsRecord output;
  output.start_idx = start_idx;
  output.batch_num = batch_num;
  output.seqs = seqs;

  return output;
#endif
}
#endif

void SeqsRead::compute() {

  int num_seqs_produced = 0;

  while (true) {
    uint64_t start_ts = getUs();

    // Read from file input, get mem_chains
    int batch_num = 0;
    bseq1_t *seqs = bseq_read(10000000, 
        &batch_num, aux->ks, aux->ks2);

    if (!seqs) break;

		if (!aux->copy_comment) {
			for (int i = 0; i < batch_num; i++) {
				free(seqs[i].comment);
				seqs[i].comment = 0;
			}
      DLOG_IF(INFO, FLAGS_v >= 2 && num_seqs_produced == 0) << "Do not append seq comment";
    }

    DLOG_IF(INFO, VLOG_IS_ON(1)) << "Read " << batch_num << " seqs in "
            << getUs() - start_ts << " us";

    // Construct output record
    SeqsRecord record;
    record.start_idx = num_seqs_produced;
    record.batch_num = batch_num;
    record.seqs = seqs;

    pushOutput(record);
    num_seqs_produced += batch_num;
  }
}

SeqsRecord SeqsToSams::compute(SeqsRecord const & input) {

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started SeqsToSams() for one input";
  uint64_t start_ts = getUs();

  bseq1_t* seqs = input.seqs;
  uint64_t start_idx = input.start_idx;
  int batch_num = input.batch_num;

  mem_alnreg_v* alnreg = new mem_alnreg_v[batch_num];

  for (int i = 0; i < batch_num; i++) {
    mem_chain_v chains = seq2chain(aux, &seqs[i]);
    kv_init(alnreg[i]);
    for (int j = 0; j < chains.n; j++) {
      mem_chain2aln(
          aux->opt, 
          aux->idx->bns, 
          aux->idx->pac,
          seqs[i].l_seq,
          (uint8_t*)seqs[i].seq,
          &chains.a[j],
          &alnreg[i]);

      free(chains.a[j].seeds);
    }
    free(chains.a);

    // Post-process each chain before output
    alnreg[i].n = mem_sort_dedup_patch(
        aux->opt, 
        aux->idx->bns, 
        aux->idx->pac, 
        (uint8_t*)seqs[i].seq, 
        alnreg[i].n, 
        alnreg[i].a);

    for (int j = 0; j < alnreg[i].n; j++) {
      mem_alnreg_t *p = &alnreg[i].a[j];
      if (p->rid >= 0 && aux->idx->bns->anns[p->rid].is_alt)
        p->is_alt = 1;
    }
  }
  mem_pestat_t pes[4];
  mem_pestat(aux->opt, aux->idx->bns->l_pac, batch_num, alnreg, pes);

#ifdef USE_HTSLIB
  for (int i =0; i< batch_num/2; i++) {
    seqs[i<<1].bams = bams_init();
    seqs[1+(i<<1)].bams = bams_init();
    mem_sam_pe(
        aux->opt,
        aux->idx->bns,
        aux->idx->pac,
        pes,
        (start_idx>>1)+i,
        &seqs[i<<1],
        &alnreg[i<<1],
        aux->h);
     }
#else
  for (int i = 0; i < batch_num/2; i++) {
    mem_sam_pe(
        aux->opt,
        aux->idx->bns,
        aux->idx->pac,
        pes,
        (start_idx>>1)+i,
        &seqs[i<<1],
        &alnreg[i<<1]);
  }
#endif
  freeAligns(alnreg, batch_num);

  // Free fields in seq except sam
  for (int i = 0; i < batch_num; i++) {
    free(seqs[i].name);
    free(seqs[i].comment);
    free(seqs[i].seq);
    free(seqs[i].qual);
  }

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished SeqsToSams() for one batch in "
    << getUs() - start_ts << " us";

  return input;
}

ChainsRecord SeqsToChains::compute(SeqsRecord const & seqs_record) {

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started SeqsToChains() for one input";

  uint64_t start_ts = getUs();
  uint64_t ref_time = 0;

  bseq1_t* seqs = seqs_record.seqs;
  uint64_t start_idx = seqs_record.start_idx;
  int batch_num = seqs_record.batch_num;

  mem_chain_v* chains = (mem_chain_v*)malloc(batch_num*sizeof(mem_chain_v));
  mem_alnreg_v* alnreg = (mem_alnreg_v*)malloc(batch_num*sizeof(mem_alnreg_v));

  for (int i = 0; i < batch_num; i++) {
    chains[i] = seq2chain(aux, &seqs[i]);
    kv_init(alnreg[i]);
  }

  ChainsRecord ret;
  ret.start_idx    = seqs_record.start_idx;
  ret.batch_num    = batch_num;
  ret.seqs         = seqs;
  ret.chains       = chains;
  ret.alnreg       = alnreg;
 
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished SeqsToChains() in " 
    << getUs() - start_ts << " us";

  return ret;
}

RegionsRecord ChainsToRegions::compute(ChainsRecord const & record) {

  uint64_t start_ts = getUs();
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started ChainsToRegions() on CPU";

  bseq1_t* seqs       = record.seqs;
  mem_chain_v* chains = record.chains;
  int batch_num       = record.batch_num;

  mem_alnreg_v* alnreg = record.alnreg;

  for (int i = 0; i < batch_num; i++) {
    for (int j = 0; j < chains[i].n; j++) {
      mem_chain2aln(
          aux->opt, 
          aux->idx->bns, 
          aux->idx->pac,
          seqs[i].l_seq,
          (uint8_t*)seqs[i].seq,
          &chains[i].a[j],
          alnreg+i);
    }
  }

  RegionsRecord output;
  output.start_idx = record.start_idx;
  output.batch_num = batch_num;
  output.seqs = seqs;
  output.alnreg = alnreg;
  output.chains = NULL;

  freeChains(chains, batch_num);

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished ChainsToRegions() on CPU for "
    << getUs() - start_ts << " us";

  return output;
}

SeqsRecord RegionsToSam::compute(RegionsRecord const & record) {

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started RegionsToSam() for one input";

  uint64_t start_ts = getUs();
  uint64_t seedcov_time = 0;

  uint64_t start_idx        = record.start_idx;
  int batch_num        = record.batch_num;
  mem_chain_v* chains  = record.chains;
  mem_alnreg_v* alnreg = record.alnreg;
  bseq1_t* seqs        = record.seqs;

  for (int i = 0; i < batch_num; i++) {
    // Post-process each chain before output
    alnreg[i].n = mem_sort_dedup_patch(
        aux->opt, 
        aux->idx->bns, 
        aux->idx->pac, 
        (uint8_t*)seqs[i].seq, 
        alnreg[i].n, 
        alnreg[i].a);

    for (int j = 0; j < alnreg[i].n; j++) {
      mem_alnreg_t *p = &alnreg[i].a[j];
      if (p->rid >= 0 && aux->idx->bns->anns[p->rid].is_alt)
        p->is_alt = 1;
    }
  }

  mem_pestat_t pes[4];
  mem_pestat(aux->opt, aux->idx->bns->l_pac, batch_num, alnreg, pes);
#ifdef USE_HTSLIB
  for (int i =0; i< batch_num/2; i++) {
    seqs[i<<1].bams = bams_init();
    seqs[1+(i<<1)].bams = bams_init();
    mem_sam_pe(
        aux->opt,
        aux->idx->bns,
        aux->idx->pac,
        pes,
        (start_idx>>1)+i,
        &seqs[i<<1],
        &alnreg[i<<1],
        aux->h);
     }
#else
  for (int i = 0; i < batch_num/2; i++) {
    mem_sam_pe(
        aux->opt,
        aux->idx->bns,
        aux->idx->pac,
        pes,
        (start_idx>>1)+i,
        &seqs[i<<1],
        &alnreg[i<<1]);
  }
#endif

  freeAligns(alnreg, batch_num);

  // Free fields in seq except sam
  for (int i = 0; i < batch_num; i++) {
    free(seqs[i].name);
    free(seqs[i].comment);
    free(seqs[i].seq);
    free(seqs[i].qual);
  }

  SeqsRecord output;
  output.start_idx = start_idx;
  output.batch_num = batch_num;
  output.seqs = seqs;

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished RegionsToSam() in " << getUs() - start_ts << " us";
  return output;
}

typedef bam1_t *bam1_p;
KSORT_INIT(sort, bam1_p, bam1_lt)

void SamsReorder::compute(int wid) {

  uint64_t n_processed = 0;
  std::unordered_map<uint64_t, SeqsRecord> record_buf;
#ifdef USE_HTSLIB
  int max_batch_records = FLAGS_max_batch_records;
  bam1_t** bam_buffer = NULL;
  int bam_buffer_idx = 0;
  int batch_records = 0;
  int bam_buffer_order = 0;
  BamsRecord output;
#else
  SeqsRecord output;
#endif
  while (true) {
    SeqsRecord input;
    SeqsRecord record;
    bool ready = this->getInput(input);
    while (!this->isFinal() && !ready) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(100));
      ready = this->getInput(input);
    }
    if (!ready) { 
      // this means isFinal() is true and input queue is empty
      break; 
    }
    if (FLAGS_inorder_output) {
      record_buf[input.start_idx] = input;

      while (record_buf.count(n_processed)) {
        record = record_buf[n_processed];
        record_buf.erase(n_processed);
        n_processed += record.batch_num;
#ifdef USE_HTSLIB
        int batch_num = record.batch_num;
        bseq1_t* seqs = record.seqs;
    
        if(!bam_buffer) {
          bam_buffer = (bam1_t**)malloc((max_batch_records*batch_num*2 )
              *sizeof(bam1_t*));
        }
        for(int i = 0; i < batch_num; i++) {
          if (seqs[i].bams) {
            for (int j =0; j < seqs[i].bams->l; j++) {
              bam_buffer[bam_buffer_idx++] = seqs[i].bams->bams[j];
            }
          }
          free(seqs[i].bams->bams);
          free(seqs[i].bams); seqs[i].bams = NULL;    
        }
        free(seqs);
 
        if (++batch_records >= max_batch_records) {
          if(FLAGS_sort) {
            uint64_t start_ts_st = getUs();
            //std::sort(bam_buffer, bam_buffer+bam_buffer_idx, bam1_lt);
            ks_mergesort(sort, bam_buffer_idx, (bam1_p *)bam_buffer, 0);
            DLOG_IF(INFO, VLOG_IS_ON(3)) << "Sorted " << bam_buffer_idx << " records in " <<
              getUs() - start_ts_st << " us";
          }
          output.bam_buffer = bam_buffer;
          output.bam_buffer_idx = bam_buffer_idx;
          output.bam_buffer_order = bam_buffer_order;
          bam_buffer_idx = 0;
          batch_records = 0;
          bam_buffer_order = bam_buffer_order + 1;
          pushOutput(output);
          bam_buffer = NULL;
        }
#else
        output = record;
        pushOutput(output);
#endif
      }
    }
    else {
      record = input;
#ifdef USE_HTSLIB
      int batch_num = record.batch_num;
      bseq1_t* seqs = record.seqs;
  
      if(!bam_buffer) {
        // in case bams->l is not 1 
        bam_buffer = (bam1_t**)malloc((max_batch_records*batch_num*2 )
            *sizeof(bam1_t*));
      }
      for(int i = 0; i < batch_num; i++) {
        if (seqs[i].bams) {
          for (int j =0; j < seqs[i].bams->l; j++) {
            bam_buffer[bam_buffer_idx++] = seqs[i].bams->bams[j];
          }
        }
        free(seqs[i].bams->bams);
        free(seqs[i].bams); seqs[i].bams = NULL;    
      }
      free(seqs);
      if (++ batch_records >=  max_batch_records) {
        if(FLAGS_sort) {
          uint64_t start_ts_st = getUs();
          //std::sort(bam_buffer, bam_buffer+bam_buffer_idx, bam1_lt);
          ks_mergesort(sort, bam_buffer_idx, (bam1_p *)bam_buffer, 0);
          DLOG_IF(INFO, VLOG_IS_ON(3)) << "Sorted " << bam_buffer_idx << " records in " <<
            getUs() - start_ts_st << " us";
        }
        output.bam_buffer = bam_buffer;
        output.bam_buffer_idx = bam_buffer_idx;
        output.bam_buffer_order = bam_buffer_order;
        bam_buffer_idx = 0;
        batch_records = 0;
        bam_buffer_order = bam_buffer_order + 1;
        pushOutput(output);
        bam_buffer = NULL;
      }
#else
      output = record;
      pushOutput(output);
#endif
    }
  }
#ifdef USE_HTSLIB
  //finish the remaining sam
  if(batch_records > 0){
    if(FLAGS_sort) {
      uint64_t start_ts_st = getUs();
      //std::sort(bam_buffer, bam_buffer+bam_buffer_idx, bam1_lt);
      ks_mergesort(sort, bam_buffer_idx, (bam1_p *)bam_buffer, 0);
      DLOG_IF(INFO, VLOG_IS_ON(3)) << "Sorted " << bam_buffer_idx << " records in " <<
        getUs() - start_ts_st << " us";
    }
    output.bam_buffer = bam_buffer;
    output.bam_buffer_idx = bam_buffer_idx;
    output.bam_buffer_order = bam_buffer_order;
    bam_buffer_idx = 0;
    batch_records = 0;
    bam_buffer_order = bam_buffer_order + 1;
    pushOutput(output);
    bam_buffer = NULL; 
  }
#endif
}

#ifdef USE_HTSLIB
int WriteOutput::compute(BamsRecord const &input)
#else
int WriteOutput::compute(SeqsRecord const &input)
#endif
{
  boost::any var = this->getConst("sam_dir");
  std::string out_dir = boost::any_cast<std::string>(var);
  bool use_file = !out_dir.empty();
#ifdef USE_HTSLIB
  // write bam output
  if (!use_file) {
    LOG(ERROR) << "Bams output only works with file output,"
      << " please specify --output_dir";
    exit(1);
  }
  const char *modes[] = {"wb", "wb0", "w"};
  samFile *fout = NULL;

  uint64_t start_ts = getUs();
  int file_id = input.bam_buffer_order;
  int bam_buffer_idx = input.bam_buffer_idx;
  bam1_t** bam_buffer = input.bam_buffer;

  std::stringstream ss;
  ss << out_dir << "/part-" << std::setw(6) << std::setfill('0') << file_id;
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "Writting to " << ss.str();
  fout = sam_open(ss.str().c_str(), modes[FLAGS_output_flag]); 
  if (!fout) {
    throw std::runtime_error("Cannot open bam output file");
  }
  int status = sam_hdr_write(fout, aux->h);
  if (status) {
    LOG(ERROR) << "sam_hdr_write error: " << status;
  }
  // start writing to file
  for (int i = 0; i < bam_buffer_idx; ++i){
    sam_write1(fout, aux->h, bam_buffer[i]); 
    bam_destroy1(bam_buffer[i]);
  }
  sam_close(fout);
  free(bam_buffer);
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Written " << bam_buffer_idx
          << " records in "
          << getUs() - start_ts << " us";
  return 0;
#else
  // write sam output
  FILE* fout;
  fout = stdout;
  uint64_t start_ts = getUs();
  int batch_num = input.batch_num;
  bseq1_t* seqs = input.seqs;
  for (int i = 0; i < batch_num; ++i) {
    if (seqs[i].sam) fputs(seqs[i].sam, fout);
    free(seqs[i].sam);
  }
  free(seqs);
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Written batch " << input.start_idx 
    << " to file in " << getUs() - start_ts << " us";
#endif
}

#ifdef USE_HTSLIB
void SamsPrint::sortAndWriteBamBatch(
    bam1_t** buf,
    int n_elements,
    std::string out_dir,
    int wid) 
{
  uint64_t start_ts = getUs();

  // sort the buffer first
  std::sort(buf, buf+n_elements, bam1_lt);

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Sort " << n_elements 
          << " records in "
          << getUs() - start_ts << " us";

  start_ts = getUs();

  bool use_file = !out_dir.empty();

  // open file if necessary
  if (use_file) {
    const char *modes[] = {"wb", "wb0", "w"};

    std::stringstream ss;
    ss << out_dir << "/part-" << wid
       << std::setw(6) << std::setfill('0') << file_id_[wid];

    fout_[wid] = sam_open(ss.str().c_str(), modes[FLAGS_output_flag]); 
    if (!fout_[wid]) {
      throw std::runtime_error("Cannot open output file");
    }
    int status = sam_hdr_write(fout_[wid], aux->h);
  }
  else {
    LOG(ERROR) << "Sorting only works with file output,"
               << " must specify --output_dir";
    exit(1);
  }

  // start writing to file
  for (int i = 0; i < n_elements; ++i){
    int status = sam_write1(fout_[wid], aux->h, buf[i]); 
    bam_destroy1(buf[i]);
  }
  if (use_file) {
    sam_close(fout_[wid]);
    file_id_[wid]++;
  }

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Written " << n_elements 
          << " records in "
          << getUs() - start_ts << " us";
}
#endif

void SamsPrint::compute(int wid) {

  boost::any var = this->getConst("sam_dir");
  std::string out_dir = boost::any_cast<std::string>(var);

  // Parameters for Sam output
  const int num_batch_per_batch = 8;
  int file_counter = 0;
  int file_id = 0;
  // If out_dir is not defined, output to stdout
  bool use_file = !out_dir.empty();

  uint64_t n_processed = 0;
  std::unordered_map<uint64_t, SeqsRecord> record_buf;

  // Open first file if output is file
#ifdef USE_HTSLIB
  samFile *fout = NULL;
  // TODO(mhhuang): To be fixed. "wbu" results in bam file header read error.
  const char *modes[] = {"wb", "wb0", "w"};
#else
  FILE* fout;
#endif

  if (!FLAGS_sort) {
    if (use_file) {
      std::stringstream ss;
      ss << out_dir  << "/part-" << wid 
        << std::setw(6) << std::setfill('0') << file_id;
#ifdef USE_HTSLIB
      DLOG_IF(INFO, VLOG_IS_ON(2)) << "Writting to " << ss.str();
      fout = sam_open(ss.str().c_str(), modes[FLAGS_output_flag]); 
      if (!fout) {
        throw std::runtime_error("Cannot open bam output file");
      }
      int status = sam_hdr_write(fout, aux->h);
#else
      fout = fopen(ss.str().c_str(), "w+");
#endif
      if (!fout) {
        throw std::runtime_error("Cannot open sam output file");
      }

      DLOG(INFO) << "Start writing output to " << ss.str();
    }
    else {
#ifdef USE_HTSLIB
      fout  = sam_open("-", modes[FLAGS_output_flag]); 
      if (!fout) {
        throw std::runtime_error("Cannot open sam output file");
      }
      int status = sam_hdr_write(fout, aux->h);
      if (status) {
        LOG(ERROR) << "sam_hdr_write error: " << status;
      }
#else
      fout = stdout;
#endif
    }
  }

  // Buffer to sort output bam
  int      max_bam_records = FLAGS_max_num_records;
#ifdef USE_HTSLIB
  bam1_t** bam_buffer;
  int      bam_buffer_idx = 0;
  if (FLAGS_sort) {
    bam_buffer = (bam1_t**)malloc(FLAGS_max_num_records*sizeof(bam1_t*));
  }
#endif

  // NOTE: input may be out-of-order, use a reorder buffer if
  // the output needs to be in-order
  while (true) {
    SeqsRecord input;
    bool ready = this->getInput(input);

    while (!this->isFinal() && !ready) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(100));
      ready = this->getInput(input);
    }
    if (!ready) { 
      // this means isFinal() is true and input queue is empty
      break; 
    }

    if (FLAGS_inorder_output) {
      // Add the current input to buffer first
      record_buf[input.start_idx] = input;

      // Find the next batch in the buffer
      while (record_buf.count(n_processed)) {
        uint64_t start_ts = getUs();

        SeqsRecord record = record_buf[n_processed];

        int      batch_num = record.batch_num;
        bseq1_t* seqs      = record.seqs;

#ifdef USE_HTSLIB    
        for (int i = 0; i < batch_num; ++i){
          if (seqs[i].bams) {   
            for (int j = 0; j < seqs[i].bams->l; j++) {   
              int status = sam_write1(fout, aux->h, seqs[i].bams->bams[j]);    
              LOG_IF(ERROR, status<0) << "sam_write1 error: " << status;
            }   
          }   
          bams_destroy(seqs[i].bams); seqs[i].bams = NULL;    
        }
#else
        for (int i = 0; i < batch_num; ++i) {
          if (seqs[i].sam) fputs(seqs[i].sam, fout);
          //err_fputs(seqs[i].sam, stdout);
          free(seqs[i].sam);
        }
#endif
        // Remove the record from buffer
        record_buf.erase(n_processed);

        n_processed += batch_num;

        DLOG_IF(INFO, VLOG_IS_ON(1)) << "Written batch " << record.start_idx << " to file in "
          << getUs() - start_ts << " us";
      }
    }
    else if (FLAGS_sort) {
#ifdef USE_HTSLIB
      uint64_t start_ts = getUs();
      int batch_num = input.batch_num;
      bseq1_t* seqs = input.seqs;

      for (int i = 0; i < batch_num; i++) {
        if (seqs[i].bams) {   
          for (int j = 0; j < seqs[i].bams->l; j++) {
            bam_buffer[bam_buffer_idx++] = seqs[i].bams->bams[j];
            if (bam_buffer_idx >= max_bam_records) {
              sortAndWriteBamBatch(bam_buffer, max_bam_records, out_dir, wid);
              bam_buffer_idx = 0;
            }
          }  
          free(seqs[i].bams->bams);
          free(seqs[i].bams); seqs[i].bams = NULL;    
        }
      }
      free(seqs);
#endif
    }
    else {
      uint64_t start_ts = getUs();
      int batch_num = input.batch_num;
      bseq1_t* seqs = input.seqs;

      if (use_file) {
        if (file_counter < num_batch_per_batch) {
          file_counter++;
        }
        else {
          file_id++;
          file_counter = 0;

          // Open a new file
          std::stringstream ss;
          ss << out_dir << "/part-" << wid
            << std::setw(6) << std::setfill('0') << file_id;

#ifdef USE_HTSLIB
          sam_close(fout);
          fout = sam_open(ss.str().c_str(), modes[FLAGS_output_flag]); 
          int status = sam_hdr_write(fout, aux->h);
#else
          fclose(fout);
          fout = fopen(ss.str().c_str(), "w+");
#endif
          DLOG(INFO) << "Start writing output to " << ss.str();
        }
      }
#ifdef USE_HTSLIB    
      for (int i = 0; i < batch_num; ++i){
        if (seqs[i].bams) {   
          for (int j = 0; j < seqs[i].bams->l; j++) {   
            int status = sam_write1(fout, aux->h, seqs[i].bams->bams[j]);    
            LOG_IF(ERROR, status<0) << "sam_write1 error: " << status;
          }   
        }   
        bams_destroy(seqs[i].bams); seqs[i].bams = NULL;    
      }
#else
      for (int i = 0; i < batch_num; ++i) {
        if (seqs[i].sam) fputs(seqs[i].sam, fout);
        //err_fputs(seqs[i].sam, stdout);
        free(seqs[i].sam);
      }
#endif
      free(seqs);

      DLOG_IF(INFO, VLOG_IS_ON(1)) << "Written batch " << input.start_idx 
        << " to file in " << getUs() - start_ts << " us";
    }
  }
#ifdef USE_HTSLIB    
  if (bam_buffer_idx > 0) {
    sortAndWriteBamBatch(bam_buffer, bam_buffer_idx, out_dir, wid);
    free(bam_buffer);
  }
  if(!FLAGS_sort) {
    sam_close(fout);
  }
  //bam_hdr_destroy(aux->h);
#else
  if (use_file) {
    fclose(fout);
  }
#endif
}
