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

#include "config.h"
#include "Pipeline.h"  
#include "util.h"
#include "bwa_wrapper.h"

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

      DLOG_IF(INFO, FLAGS_v >= 2) << "Serializing seq batch of "
        << length / 1024  << "kb"
        << " in " << getUs() - start_ts << " us";

      start_ts = getUs();

      // First query a process to send data to
      int proc_id = -1;
      bwaMPIRecv(&proc_id, 1, MPI::INT, MPI::ANY_SOURCE, SEQ_DP_QUERY);

      bwaMPISend(&length, 1, MPI::INT, proc_id, SEQ_DP_LENGTH);

      bwaMPISend(ser_data.c_str(), length, MPI::CHAR, proc_id, SEQ_DP_DATA);

      DLOG_IF(INFO, FLAGS_v >= 1) << "Sending seqs batch " << record.start_idx
        << " to proc_" << proc_id
        << " takes " << getUs() - start_ts << " us";
    }
    else {
      // this means isFinal() is true and input queue is empty
      DLOG_IF(INFO, FLAGS_v >= 1) << "Finish reading seqs, start send finish signals";

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

      DLOG_IF(INFO, FLAGS_v >= 1) << "Receive one read batch in "
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
      DLOG_IF(INFO, FLAGS_v >= 2) << "Serializing sam batch of "
        << length / 1024  << "kb"
        << " in " << getUs() - start_ts << " us";

      start_ts = getUs();

      // Send proc_id to master to let master receive following msg
      bwaMPISend(&rank, 1, MPI::INT, MASTER_RANK, SAM_RV_QUERY);

      bwaMPISend(&length, 1, MPI::INT, MASTER_RANK, SAM_RV_LENGTH);

      bwaMPISend(ser_data.c_str(), length,
          MPI::CHAR, MASTER_RANK, SAM_RV_DATA);

      DLOG_IF(INFO, FLAGS_v >= 1) << "Sending sam batch " << input.start_idx
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

    DLOG_IF(INFO, FLAGS_v >= 1) << "Read " << batch_num << " seqs in "
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

  DLOG_IF(INFO, FLAGS_v >= 1) << "Started SeqsToSams() for one input";
  uint64_t start_ts = getUs();

  bseq1_t* seqs = input.seqs;
  int start_idx = input.start_idx;
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

  DLOG_IF(INFO, FLAGS_v >= 1) << "Finished SeqsToSams() for one batch in "
    << getUs() - start_ts << " us";

  return input;
}

ChainsRecord SeqsToChains::compute(SeqsRecord const & seqs_record) {

  DLOG_IF(INFO, FLAGS_v >= 1) << "Started SeqsToChains() for one input";

  uint64_t start_ts = getUs();
  uint64_t ref_time = 0;

  bseq1_t* seqs = seqs_record.seqs;
  int start_idx = seqs_record.start_idx;
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
 
  DLOG_IF(INFO, FLAGS_v >= 1) << "Finished SeqsToChains() in " << getUs() - start_ts << " us";

  return ret;
}

inline int cal_max_gap(
    const mem_opt_t *opt, 
    int qlen
) {
	int l_del = (int)((double)(qlen * opt->a - opt->o_del) / opt->e_del + 1.);
	int l_ins = (int)((double)(qlen * opt->a - opt->o_ins) / opt->e_ins + 1.);
	int l = l_del > l_ins? l_del : l_ins;
	l = l > 1? l : 1;
	return l < opt->w<<1? l : opt->w<<1;
}

inline void getChainRef(
      const ktp_aux_t* aux,
      const bseq1_t* seq,
      const mem_chain_v* chain,
      mem_chainref_t* &ref
) {
  int64_t l_pac = aux->idx->bns->l_pac;

  // input for SmithWaterman for each chain
  ref = (mem_chainref_t*)malloc(chain->n*sizeof(mem_chainref_t));

  for (int j = 0; j < chain->n; j++) {
    // Prepare the maxspan and rseq for each seed
    mem_chain_t *c = &chain->a[j]; 

    int64_t rmax[2], tmp, max = 0;
    uint8_t *rseq = 0;
    uint64_t *srt;

    if (c->n == 0) {
      continue;
    }
    // get the max possible span
    rmax[0] = l_pac<<1; rmax[1] = 0;
    for (int k = 0; k < c->n; ++k) {
      int64_t b, e;
      const mem_seed_t *t = &c->seeds[k];
      b = t->rbeg - (t->qbeg + cal_max_gap(aux->opt, t->qbeg));
      e = t->rbeg + t->len + ((seq->l_seq - t->qbeg - t->len)
          + cal_max_gap(aux->opt, seq->l_seq - t->qbeg - t->len));
      rmax[0] = rmax[0] < b? rmax[0] : b;
      rmax[1] = rmax[1] > e? rmax[1] : e;
      if (t->len > max) max = t->len;
    }
    rmax[0] = rmax[0] > 0? rmax[0] : 0;
    rmax[1] = rmax[1] < l_pac<<1? rmax[1] : l_pac<<1;
    if (rmax[0] < l_pac && l_pac < rmax[1]) { // crossing the forward-reverse boundary; then choose one side
      if (c->seeds[0].rbeg < l_pac) rmax[1] = l_pac; // this works because all seeds are guaranteed to be on the same strand
      else rmax[0] = l_pac;
    }
    // retrieve the reference sequence
    int rid;
    rseq = bns_fetch_seq(aux->idx->bns, aux->idx->pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);
    assert(c->rid == rid);

    srt = (uint64_t *)malloc(c->n * 8);
    for (int l = 0; l < c->n; ++l)
      srt[l] = (uint64_t)c->seeds[l].score<<32 | l;
    ks_introsort_64(c->n, srt);

    ref[j].rmax[0] = rmax[0];
    ref[j].rmax[1] = rmax[1];
    ref[j].rseq = rseq;
    ref[j].srt = srt;
  }
}

inline void packReadData(ktp_aux_t* aux, 
    const bseq1_t* seq, 
    const mem_chain_v* chains,
    mem_alnreg_v* alnregs, 
    char* buffer, 
    int &buffer_idx, 
    int &task_num, 
    mem_alnreg_t** &region_batch,
    mem_chain_t** &chain_batch) 
{
  mem_chainref_t* ref;
  int chain_idx = 0;
  int seed_idx = chains->a[0].n - 1;
  getChainRef(aux, seq, chains, ref);

  for ( ; chain_idx < chains->n; 
          seed_idx > 0 ? seed_idx-- : seed_idx = chains->a[++chain_idx].n-1)
  {
      uint32_t sorted_idx = (uint32_t)(ref[chain_idx].srt[seed_idx]);

    // get next available seed in the current read
      mem_seed_t* seed_array = &chains->a[chain_idx].seeds[sorted_idx];
      // initialize the newreg
      mem_alnreg_t* newreg = kv_pushp(mem_alnreg_t, *alnregs);
      memset(newreg, 0, sizeof(mem_alnreg_t));

      newreg->score  = seed_array->len * aux->opt->a;
      newreg->truesc = seed_array->len * aux->opt->a;
      newreg->qb = seed_array->qbeg;
      newreg->rb = seed_array->rbeg;
      newreg->qe = seed_array->qbeg + seed_array->len;
      newreg->re = seed_array->rbeg + seed_array->len; 
      newreg->rid = chains->a[chain_idx].rid;
      newreg->seedlen0 = seed_array->len; 
      newreg->frac_rep = chains->a[chain_idx].frac_rep;
      newreg->w = aux->opt->w;
  }
  chain_idx = 0;
  seed_idx = 0;
  int counter8 = 0;
  int tmp_int = 0;
  int chain_num = 0;
  int seed_num = 0;
  int chain_num_addr = 0;
  int seed_num_addr = 0;
  int idx_end_addr = buffer_idx;
  buffer_idx += 4;

  // pack the read sequence
  *((int*)(&buffer[buffer_idx])) = seq->l_seq;
  buffer_idx += 4;
  for ( int i = 0; i < seq->l_seq; i++ ){
    counter8 += 1;
    tmp_int = tmp_int << 4 | seq->seq[i];
    if ( counter8 % 8 ==0 ) {
      *(uint32_t*) (&buffer[buffer_idx]) = tmp_int ;
      buffer_idx += 4 ;
    }
  }
  if ( counter8 % 8 !=0 ) {
    tmp_int = tmp_int << (4*(8 - counter8 % 8));
    *(uint32_t*) (&buffer[buffer_idx]) = tmp_int ;
    buffer_idx += 4 ;
  }
  counter8 = 0;
  tmp_int = 0;
  chain_num_addr = buffer_idx;
  buffer_idx += 4;
  int region_num = 0; // used to keep track of the current region

  for ( chain_idx = 0; chain_idx < chains->n; chain_idx++) {
    // Pack the maxspan and rseq
    *((int64_t*)(&buffer[buffer_idx]))= ref[chain_idx].rmax[0];
    buffer_idx += 8;
    *((int64_t*)(&buffer[buffer_idx]))= ref[chain_idx].rmax[1];
    buffer_idx += 8;
    for ( int i = 0; i < ref[chain_idx].rmax[1] - ref[chain_idx].rmax[0]; i++) {
      counter8 = counter8 + 1;
      tmp_int = tmp_int << 4 | ref[chain_idx].rseq[i];
      if ( counter8 % 8 ==0 ) {
        *(uint32_t*) (&buffer[buffer_idx]) = tmp_int ;
        buffer_idx += 4 ;
      }
    } 
    if ( counter8 % 8 !=0 ) {
      tmp_int = tmp_int << (4*(8 - counter8 % 8));
      *(uint32_t*) (&buffer[buffer_idx]) = tmp_int ;
      buffer_idx += 4 ;
    }
    counter8 = 0;
    // record the address of seed number
    seed_num_addr = buffer_idx ;
    seed_num = 0;
    buffer_idx += 4 ; 
    // pack the seed information
    for ( seed_idx = chains->a[chain_idx].n -1 ; seed_idx >= 0 ; seed_idx--) {
      uint32_t sorted_idx = (uint32_t)(ref[chain_idx].srt[seed_idx]);
      // get next available seed in the current read
      mem_seed_t* seed_array = &chains->a[chain_idx].seeds[sorted_idx];
     
      if (seed_array->qbeg > 0 ||
          seed_array->qbeg + seed_array->len != seq->l_seq)
      {
        region_batch[task_num] = &(alnregs->a[region_num]) ;
        chain_batch[task_num] = &(chains->a[chain_idx]) ;
        *((int*)(&buffer[buffer_idx])) = task_num ; 
        buffer_idx += 4;
        task_num += 1;
        seed_num += 1;
        *((int64_t*)(&buffer[buffer_idx])) = seed_array->rbeg ;
        buffer_idx += 8;
        *((int32_t*)(&buffer[buffer_idx])) = seed_array->qbeg ;
        buffer_idx += 4;
        *((int32_t*)(&buffer[buffer_idx])) = seed_array->len ;
        buffer_idx += 4;
      }
      else {
        // no need to loop again 
        alnregs->a[region_num].seedcov = seed_array->len;
      }
      region_num++;
    }
     *((int*)(&buffer[seed_num_addr])) = seed_num ;
     chain_num += 1;
  }
  *((int*)(&buffer[chain_num_addr])) = chain_num ;
  *((int*)(&buffer[idx_end_addr])) = buffer_idx/4 ;

  for (int i=0; i< chains->n; i++) {
    free(ref[i].srt);
    free(ref[i].rseq);
  }
  free(ref);
}

void ChainsToRegionsFPGA::compute(int wid) {

  int chunk_size = FLAGS_chunk_size;

  const int stage_num = 2;
  int stage_cnt = 0;

  // Create FPGAAgent
  FPGAAgent agent(opencl_env, chunk_size);

  int task_num = 0;

  // Batch of regions
  mem_alnreg_t** region_batch[stage_num];
  mem_chain_t** chain_batch[stage_num];
  int stage_task_num_a[stage_num];
  int stage_task_num_b[stage_num];
  for (int i = 0; i < stage_num; i++) {
    region_batch[i] = new mem_alnreg_t*[2*chunk_size];
    chain_batch[i] = new mem_chain_t*[2*chunk_size];
    stage_task_num_a[i] = 0;
    stage_task_num_b[i] = 0;
  } 
  // elements of the Regions record
  bseq1_t* seqs        ;
  mem_chain_v* chains  ;
  mem_alnreg_v* alnreg ;

  // kernel_buffer 
  char* kernel_buffer = new char [10000*chunk_size];
  short* kernel_buffer_out = new short [40 * chunk_size];
  int kernel_buffer_idx = 0;
  // For statistics
  uint64_t start_ts;
  uint64_t last_batch_ts;
  uint64_t last_output_ts;
  uint64_t last_read_ts;

  uint64_t process_time  = 0;
  uint64_t pending_time  = 0;
  uint64_t finish_time   = 0;
  uint64_t nextTask_time = 0;
  uint64_t batch_time    = 0;
  uint64_t output_time   = 0;

  uint64_t nextTask_num = 0;
  int      start_idx    = 0;
  int      batch_num    = 0;

  bool flag_need_reads = false;
  bool flag_more_reads = true;
  ChainsRecord record; 
  while (flag_more_reads) { 
    // get one Chains record
    bool ready = this->getInput(record);
    start_ts = getUs();
    while (!this->isFinal() && !ready) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(100));
      ready = this->getInput(record);
    }
    DLOG_IF(INFO, FLAGS_v >= 2) << "Wait for input for FPGA takes " << getUs() - start_ts << " us";
    if(!ready) {
      flag_more_reads = false;
      break;
    }
    DLOG_IF(INFO, FLAGS_v >= 1) << "Started ChainsToRegions() on FPGA ";

    last_batch_ts  = getUs();
    last_output_ts = last_batch_ts;
    last_read_ts  = last_batch_ts;
    start_idx = record.start_idx;
    batch_num = record.batch_num;
    seqs = record.seqs;
    chains = record.chains;
    alnreg = record.alnreg;

    int i = 0;
    uint64_t last_prepare_ts = getUs();
    uint64_t last_batch_ts = getUs();
    int kernel_task_num_a = 0;
    int kernel_task_num_b = 0;
    int kernel_size_a = 0;
    int kernel_size_b = 0;
    bool reach_half = false;
    while (i < batch_num) {
      if (task_num < chunk_size/2) {
        packReadData(aux, seqs+i, chains+i,
            alnreg+i, kernel_buffer, kernel_buffer_idx, task_num,
            region_batch[stage_cnt], chain_batch[stage_cnt]);
        i++;
      }
      else if (task_num >= chunk_size/2 && reach_half == false) {
        kernel_size_a = kernel_buffer_idx;
        kernel_task_num_a = task_num;
        kernel_buffer_idx = 0;
        reach_half = true;
      }
      else if (task_num < chunk_size) {
        packReadData(aux, seqs+i, chains+i,
            alnreg+i, &kernel_buffer[kernel_size_a], kernel_buffer_idx, task_num,
            region_batch[stage_cnt], chain_batch[stage_cnt]);
        i++;
      }
      else if (task_num >= chunk_size) {
        kernel_size_b = kernel_buffer_idx;
        kernel_task_num_b = task_num - kernel_task_num_a;
        DLOG_IF(INFO, FLAGS_v >= 3) << "Prepare data for FPGA  "<< stage_cnt << " takes " << getUs()- last_prepare_ts << " us " ;
        stage_task_num_a[stage_cnt] = kernel_task_num_a ;
        stage_task_num_b[stage_cnt] = kernel_task_num_b ;
        //extendOnFPGA(&agent, kernel_buffer, kernel_size_a, kernel_size_b, stage_cnt);
        //stage_cnt = 1 - stage_cnt;
        //if (agent.pending(stage_cnt)){
        //  FPGAPostProcess(&agent, kernel_buffer_out, stage_task_num_a[stage_cnt], 
        //      stage_task_num_b[stage_cnt],region_batch[stage_cnt], chain_batch[stage_cnt], stage_cnt );
        //}
        stage_cnt = 1 - stage_cnt;
        if (agent.pending(stage_cnt)){
          FPGAPostProcess(&agent, kernel_buffer_out, stage_task_num_a[stage_cnt], 
              stage_task_num_b[stage_cnt],region_batch[stage_cnt], chain_batch[stage_cnt], stage_cnt );
        }
        extendOnFPGA(&agent, kernel_buffer, kernel_size_a, kernel_size_b, 1 - stage_cnt);
        last_prepare_ts = getUs();
        DLOG_IF(INFO, FLAGS_v >= 3) << "This chunk of FPGA takes " << getUs()- last_batch_ts<< " us " ;
        last_batch_ts = getUs();
        task_num = 0;
        kernel_buffer_idx = 0;
        kernel_size_a = 0;
        kernel_size_b = 0;
        kernel_task_num_a = 0;
        kernel_task_num_b = 0;
        reach_half = false;
      }
    }
    if (agent.pending(1 - stage_cnt)){
      FPGAPostProcess(&agent, kernel_buffer_out, stage_task_num_a[1 - stage_cnt], 
          stage_task_num_b[1 - stage_cnt], region_batch[1 - stage_cnt],
          chain_batch[1 - stage_cnt], 1 - stage_cnt );
    }
    // finish the remain reads even with small task number 
    extendOnFPGA(&agent, kernel_buffer, kernel_size_a, kernel_size_b, stage_cnt);
    FPGAPostProcess(&agent, kernel_buffer_out, kernel_task_num_a, kernel_task_num_b,
        region_batch[stage_cnt], chain_batch[stage_cnt], stage_cnt );
    task_num = 0;
    kernel_buffer_idx = 0;
    kernel_size_a = 0;
    kernel_size_b = 0;
    kernel_task_num_a = 0;
    kernel_task_num_b = 0;
    reach_half = false;
    i = 0;

    RegionsRecord outputRecord;
    outputRecord.start_idx = start_idx;
    outputRecord.batch_num = batch_num;
    outputRecord.seqs = seqs;
    outputRecord.chains = chains;
    outputRecord.alnreg = alnreg;

    freeChains(chains, batch_num);
    pushOutput(outputRecord);
    DLOG_IF(INFO, FLAGS_v >= 1) << "Finished ChainsToRegions() on FPGA for "
      << getUs() - last_output_ts << " us";
    last_output_ts = getUs();
  }
}

RegionsRecord ChainsToRegions::compute(ChainsRecord const & record) {

  uint64_t start_ts = getUs();
  DLOG_IF(INFO, FLAGS_v >= 1) << "Started ChainsToRegions() on CPU";

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

  DLOG_IF(INFO, FLAGS_v >= 1) << "Finished ChainsToRegions() on CPU for "
    << getUs() - start_ts << " us";

  return output;
}

SeqsRecord RegionsToSam::compute(RegionsRecord const & record) {

  DLOG_IF(INFO, FLAGS_v >= 1) << "Started RegionsToSam() for one input";

  uint64_t start_ts = getUs();
  uint64_t seedcov_time = 0;

  int start_idx        = record.start_idx;
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

  DLOG_IF(INFO, FLAGS_v >= 1) << "Finished RegionsToSam() in " << getUs() - start_ts << " us";
  return output;
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

  DLOG_IF(INFO, FLAGS_v >= 1) << "Sort " << n_elements 
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

  DLOG_IF(INFO, FLAGS_v >= 1) << "Written " << n_elements 
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
      DLOG_IF(INFO, FLAGS_v >= 2) << "Writting to " << ss.str();
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

        DLOG_IF(INFO, FLAGS_v >= 1) << "Written batch " << record.start_idx << " to file in "
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

      DLOG_IF(INFO, FLAGS_v >= 1) << "Written batch " << input.start_idx << " to file in "
        << getUs() - start_ts << " us";
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
