#include <algorithm>
#include <boost/asio.hpp>
#include <boost/function_types/result_type.hpp>
#include <boost/make_shared.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread.hpp>
#include <list>
#include <queue>

#include "bwa/utils.h"
#include "htslib/ksort.h"

#include "config.h"
#include "Pipeline.h"  
#include "util.h"
#include "bwa_wrapper.h"

// Comparator function for bam1_t records
bool bam1_lt(const bam1_t* a, const bam1_t* b) {
  return ((uint64_t)a->core.tid<<32|(a->core.pos+1)<<1|bam_is_rev(a)) 
       < ((uint64_t)b->core.tid<<32|(b->core.pos+1)<<1|bam_is_rev(b));
}

#ifdef BUILD_FPGA
#include "Extension.h"
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

      VLOG(2) << "Serializing seq batch of "
        << length / 1024  << "kb"
        << " in " << getUs() - start_ts << " us";

      start_ts = getUs();

      // First query a process to send data to
      int proc_id = -1;
      bwaMPIRecv(&proc_id, 1, MPI::INT, MPI::ANY_SOURCE, SEQ_DP_QUERY);

      bwaMPISend(&length, 1, MPI::INT, proc_id, SEQ_DP_LENGTH);

      bwaMPISend(ser_data.c_str(), length, MPI::CHAR, proc_id, SEQ_DP_DATA);

      VLOG(1) << "Sending seqs batch " << record.start_idx
        << " to proc_" << proc_id
        << " takes " << getUs() - start_ts << " us";
    }
    else {
      // this means isFinal() is true and input queue is empty
      VLOG(1) << "Finish reading seqs, start send finish signals";

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

      VLOG(1) << "Receive one read batch in "
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
      VLOG(2) << "Serializing sam batch of "
        << length / 1024  << "kb"
        << " in " << getUs() - start_ts << " us";

      start_ts = getUs();

      // Send proc_id to master to let master receive following msg
      bwaMPISend(&rank, 1, MPI::INT, MASTER_RANK, SAM_RV_QUERY);

      bwaMPISend(&length, 1, MPI::INT, MASTER_RANK, SAM_RV_LENGTH);

      bwaMPISend(ser_data.c_str(), length,
          MPI::CHAR, MASTER_RANK, SAM_RV_DATA);

      VLOG(1) << "Sending sam batch " << input.start_idx
        << " to master takes " << getUs() - start_ts << " us";

      //freeSeqs(input.seqs, input.batch_num);
      for (int i = 0; i < input.batch_num; i++) {
        free(input.seqs[i].sam);
      }
      free(input.seqs);
    }
  }
}

std::string SamsSend::serialize(SeqsRecord* data) {

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
}

void SamsReceive::compute() {
  
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
}

SeqsRecord SamsReceive::deserialize(const char* data, size_t length) {

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
}
#endif

void SeqsRead::compute() {

  int num_seqs_produced = 0;

  while (true) {
    uint64_t start_ts = getUs();

    // Read from file input, get mem_chains
    int batch_num = 0;
    bseq1_t *seqs = bseq_read(aux->actual_chunk_size, 
        &batch_num, aux->ks, aux->ks2);

    if (!seqs) break;

		if (!aux->copy_comment) {
			for (int i = 0; i < batch_num; i++) {
				free(seqs[i].comment);
				seqs[i].comment = 0;
			}
      VLOG_IF(2, num_seqs_produced == 0) << "Do not append seq comment";
    }

    VLOG(1) << "Read " << batch_num << " seqs in "
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

  VLOG(2) << "Started SeqsToSams() for one input";
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

  VLOG(1) << "Compute a batch in "
    << getUs() - start_ts << " us";

  return input;
}

inline void SeqsToChains::prepareChainRef(
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

ChainsRecord SeqsToChains::compute(SeqsRecord const & seqs_record) {

  VLOG(2) << "Started SeqsToChains() for one input";

  uint64_t start_ts = getUs();
  uint64_t ref_time = 0;

  bseq1_t* seqs = seqs_record.seqs;
  int start_idx = seqs_record.start_idx;
  int batch_num = seqs_record.batch_num;

  mem_chain_v* chains = (mem_chain_v*)malloc(batch_num*sizeof(mem_chain_v));
  mem_alnreg_v* alnreg = (mem_alnreg_v*)malloc(batch_num*sizeof(mem_alnreg_v));

  mem_chainref_t** chainrefs;
  std::list<SWRead*>* read_batch;
  std::vector<int>* chains_idxes;

#ifdef BUILD_FPGA
  if (FLAGS_use_fpga && FLAGS_max_fpga_thread) {
    // Freed one by one in chain2reg stage
    read_batch = new std::list<SWRead*>;
    // Allocate here and it will be freed after SWRead is destroyed
    chainrefs = (mem_chainref_t**)malloc(batch_num*sizeof(mem_chainref_t*));
    // Freed in reg2sam stage
    chains_idxes = new std::vector<int>[batch_num];
  }
#endif

  for (int i = 0; i < batch_num; i++) {
    chains[i] = seq2chain(aux, &seqs[i]);
    kv_init(alnreg[i]);

#ifdef BUILD_FPGA
    uint64_t start_ts = getNs();
    if (FLAGS_use_fpga && FLAGS_max_fpga_thread) {
      prepareChainRef(aux, &seqs[i], &chains[i], chainrefs[i]);

      SWRead *read_ptr = new SWRead(start_idx, i, aux, 
          seqs+i, chains+i, alnreg+i, chainrefs[i], chains_idxes+i);

      read_batch->push_back(read_ptr); 
    }
    ref_time += getNs() - start_ts;
#endif
  }

  ChainsRecord ret;
  ret.start_idx    = seqs_record.start_idx;
  ret.batch_num    = batch_num;
  ret.seqs         = seqs;
  ret.chains       = chains;
  ret.alnreg       = alnreg;
#ifdef BUILD_FPGA
  if (FLAGS_use_fpga && FLAGS_max_fpga_thread) {
    ret.read_batch = read_batch;
    ret.chains_idxes = chains_idxes;
  }
  else 
#endif
  {
    ret.read_batch   = NULL;
    ret.chains_idxes = NULL;
  }

  VLOG(1) << "Produced a chain batch in " << getUs() - start_ts << " us";
  VLOG_IF(2, FLAGS_use_fpga) << "prepareChainRef takes " << ref_time/1e3 << " us";

  return ret;
}

inline bool ChainsToRegions::addBatch(
    std::list<SWRead*> &read_batch,
    std::unordered_map<uint64_t, int> &tasks_remain,
    std::unordered_map<uint64_t, ChainsRecord> &input_buf,
    std::unordered_map<uint64_t, RegionsRecord> &output_buf) 
{
  ChainsRecord record;
  bool ready = this->getInput(record);

  uint64_t start_ts = getUs();
  while (!this->isFinal() && !ready) {
    boost::this_thread::sleep_for(boost::chrono::microseconds(100));
    ready = this->getInput(record);
  }
  VLOG(2) << "Wait for input takes " << getUs() - start_ts << " us";
  if (!ready) { 
    // this means isFinal() is true and input queue is empty
    return false; 
  }

  // Get input record
  int start_idx        = record.start_idx;
  int batch_num        = record.batch_num;
  bseq1_t* seqs        = record.seqs;
  mem_chain_v* chains  = record.chains;
  mem_alnreg_v* alnreg = record.alnreg;
  std::vector<int>* chains_idxes = record.chains_idxes;
  std::list<SWRead*>* new_reads  = record.read_batch;

  input_buf[start_idx] = record;

  // push record to output table
  RegionsRecord output;
  output.start_idx = start_idx;
  output.batch_num = batch_num;
  output.seqs = seqs;
  output.chains = chains;
  output.alnreg = alnreg;
  output.chains_idxes = chains_idxes;

  output_buf[start_idx] = output;
  tasks_remain[start_idx] = batch_num;

  // copy all new reads to current read_batch
  read_batch.splice(read_batch.end(), *new_reads);

  delete new_reads;

  return true;
}

void ChainsToRegions::compute(int wid) {

#ifdef BUILD_FPGA
  if (FLAGS_use_fpga && wid < FLAGS_max_fpga_thread) {

    VLOG(1) << "Worker " << wid << " is working on FPGA";

    int chunk_size = FLAGS_chunk_size;

    const int stage_num = 2;
    int stage_cnt = 0;

    // Create FPGAAgent
    FPGAAgent agent(opencl_env, chunk_size);

    int task_num = 0;

    // Batch of SWTasks
    ExtParam** task_batch[stage_num];

    for (int i = 0; i < stage_num; i++) {
      task_batch[i] = new ExtParam*[chunk_size];
    }

    // Batch of SWReads
    std::list<SWRead*> read_batch;

    // Table to keep track of each record
    std::unordered_map<uint64_t, int> tasks_remain;
    std::unordered_map<uint64_t, ChainsRecord> input_buf;
    std::unordered_map<uint64_t, RegionsRecord> output_buf;

    // For statistics
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
    int      batch_num    = 0;
    int      read_num     = 0;

    bool flag_need_reads = false;
    bool flag_more_reads = true;

    while (flag_more_reads || !read_batch.empty()) { 

      if (read_batch.empty()) {
        // get initial input batch
        if (!addBatch(read_batch, tasks_remain, input_buf, output_buf)) {
          flag_more_reads = false;
        }

        last_batch_ts  = getUs();
        last_output_ts = last_batch_ts;
        last_read_ts  = last_batch_ts;
      }
      else {
        std::list<SWRead*>::iterator iter = read_batch.begin();

        while (iter != read_batch.end()) {

          uint64_t start_ts;
          uint64_t wait_ts;
          uint64_t start_idx;

          ExtParam* param_task;
          start_ts = getNs();
          SWRead::TaskStatus status = (*iter)->nextTask(param_task);
          nextTask_time += getNs() - start_ts; 
          nextTask_num ++;

          int curr_size = 0;

          switch (status) {

            case SWRead::TaskStatus::Successful:

              task_batch[stage_cnt][task_num] = param_task;
              task_num++;
              if (task_num >= chunk_size) {
                VLOG(3) << "nextTask takes " << nextTask_time/1e3 << " us in " 
                        << nextTask_num << " calls";

                nextTask_time = 0;
                nextTask_num = 0;

                extendOnFPGAPackInput(
                    &agent,
                    stage_cnt,
                    task_batch[stage_cnt],
                    chunk_size,
                    aux->opt);

                stage_cnt = 1 - stage_cnt;

                // Wait for last batch to finish if valid
                if (agent.pending(stage_cnt)) {
                  extendOnFPGAProcessOutput(
                      &agent,
                      stage_cnt,
                      task_batch[stage_cnt],
                      chunk_size,
                      aux->opt);
                }

                VLOG(3) << "Process SWRead::Pending takes " << pending_time/1e3 << " us";
                VLOG(3) << "Process SWRead::Finish takes " << finish_time/1e3 << " us";
                VLOG(3) << "Batch takes " << getUs() - last_batch_ts << " us";

                pending_time = 0;
                finish_time  = 0;

                if (batch_num < 100) {
                  batch_num ++;
                  batch_time += getUs() - last_batch_ts;
                }
                else {
                  VLOG(2) << "Extension task time is "
                    << (double)batch_time / batch_num
                    << " us/chunk";
                  batch_num = 0;
                  batch_time = 0;
                }
                last_batch_ts = getUs();

                task_num = 0;
              }
              iter ++;
              break;

            case SWRead::TaskStatus::Pending:
              start_ts = getNs();
              if (flag_more_reads) {
                // Try to get a new batch
                curr_size = read_batch.size(); 

                if (addBatch(read_batch, tasks_remain, input_buf, output_buf)) {
                  iter = read_batch.begin();
                  std::advance(iter, curr_size);
                }
                else {
                  flag_more_reads = false;
                }
              }
              else {
                // No more new tasks, must do extend before proceeding
                if (agent.pending(stage_cnt)) {
                  extendOnFPGAProcessOutput(
                      &agent,
                      stage_cnt,
                      task_batch[stage_cnt],
                      chunk_size,
                      aux->opt);
                }
                else {
                  extendOnCPU(task_batch[stage_cnt], task_num, aux->opt);
                  stage_cnt = 1 - stage_cnt;

                  task_num = 0;
                  iter++;
                }
              }
              pending_time += getNs() - start_ts;
              break;

            case SWRead::TaskStatus::Finished:
              start_ts = getNs();
              // Read is finished, remove from batch
              start_idx = (*iter)->start_idx();
              tasks_remain[start_idx]--;

              // Free SWRead
              delete *iter;

              iter = read_batch.erase(iter);

              // Collecting throughput for reads
              if (read_num < 10000) {
                read_num ++;
              }
              else {
                VLOG(2) << "Finished read throughput is "
                  << (double)(getUs() - last_read_ts)/read_num
                  << " us/read";
                read_num = 0;
                last_read_ts = getUs();
              }

              // Check if corresponding batch is finished
              if (tasks_remain[start_idx] == 0) {
                pushOutput(output_buf[start_idx]);


                tasks_remain.erase(start_idx);
                input_buf.erase(start_idx);
                output_buf.erase(start_idx);

                VLOG(1) << "Produced a region batch in "
                  << getUs() - last_output_ts << " us";
                VLOG(1) << "Currently there are " << tasks_remain.size()
                  << " active batches";

                last_output_ts = getUs();
              }
              finish_time += getNs() - start_ts;
              break;

            default: ;
          }
        }
      }
    }

    for (int i = 0; i < stage_num; i++) {
      delete [] task_batch[i];
    }
  }
  else 
#endif
  {
    VLOG(1) << "Worker " << wid << " is working on CPU";

    uint64_t last_output_ts;
    uint64_t last_read_ts;

    ChainsRecord record;
    bool ready = this->getInput(record);

    uint64_t start_ts = getUs();
    if (!ready) { 
      // this means isFinal() is true or timeout
      return; 
    }
    VLOG(2) << "Wait for input takes " << getUs() - start_ts << " us";

    last_output_ts = getUs();
    last_read_ts = getUs();

    bseq1_t* seqs       = record.seqs;
    mem_chain_v* chains = record.chains;
    int batch_num       = record.batch_num;

#ifdef BUILD_FPGA
    // Free all read batches
    if (FLAGS_use_fpga && FLAGS_max_fpga_thread) {
      std::list<SWRead*>* read_batch = record.read_batch;
      for (std::list<SWRead*>::iterator iter = read_batch->begin();
          iter != read_batch->end();
          iter ++)
      {
        delete *iter;
      }
      delete read_batch;
    }
#endif

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
    VLOG(2) << "Finished read throughput is "
      << (double)(getUs() - last_read_ts)/batch_num << " us/read";

    RegionsRecord output;
    output.start_idx = record.start_idx;
    output.batch_num = batch_num;
    output.seqs = seqs;
    output.alnreg = alnreg;
    output.chains = NULL;
    output.chains_idxes = NULL;

    freeChains(chains, batch_num);
    delete [] record.chains_idxes;

    VLOG(1) << "Produced a region batch in "
      << getUs() - last_output_ts << " us";

    pushOutput(output);
  }
}

SeqsRecord RegionsToSam::compute(RegionsRecord const & record) {

  VLOG(2) << "Started RegionsToSam() for one input";

  uint64_t start_ts = getUs();
  uint64_t seedcov_time = 0;

  int start_idx        = record.start_idx;
  int batch_num        = record.batch_num;
  mem_chain_v* chains  = record.chains;
  mem_alnreg_v* alnreg = record.alnreg;
  bseq1_t* seqs        = record.seqs;
  std::vector<int>* chains_idxes = record.chains_idxes;

  for (int i = 0; i < batch_num; i++) {
#ifdef BUILD_FPGA
    uint64_t start_ts = getNs();
    if (chains_idxes) {
      // Calculate seed coverage
      for (int j = 0; j < alnreg[i].n; j++) {
        mem_alnreg_t *newreg = &alnreg[i].a[j];
        int chain_idx = chains_idxes[i][j];

        int seedcov = 0;
        for (int k = 0; k < chains[i].a[chain_idx].n; k++) {
          const mem_seed_t *seed = &chains[i].a[chain_idx].seeds[k];
          if (seed->qbeg >= newreg->qb && 
              seed->qbeg + seed->len <= newreg->qe && 
              seed->rbeg >= newreg->rb && 
              seed->rbeg + seed->len <= newreg->re){
            seedcov += seed->len; 
          }
        }
        newreg->seedcov = seedcov;  
      }
    }
    seedcov_time += getNs() - start_ts;
#endif
    
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
#ifdef BUILD_FPGA
  if (chains_idxes) {
    delete [] chains_idxes;
    freeChains(chains, batch_num);
  }
  VLOG(2) << "Seed coverage time is " << seedcov_time/1e3 << " us";
#endif

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

  VLOG(1) << "Produced a sam batch in "
    << getUs() - start_ts << " us";

  return output;
}

void SamsPrint::sortAndWriteBamBatch(
    bam1_t** buf,
    int n_elements,
    std::string out_dir) 
{
  uint64_t start_ts = getUs();

  // sort the buffer first
  std::sort(buf, buf+n_elements, bam1_lt);

  VLOG(1) << "Sort " << n_elements 
          << " records in "
          << getUs() - start_ts << " us";

  start_ts = getUs();

  bool use_file = !out_dir.empty();
  // open file if necessary
  if (use_file) {
    const char *modes[] = {"wb", "wbu", "w"};

    std::stringstream ss;
    ss << out_dir << "/part-"
       << std::setw(6) << std::setfill('0') << file_id_;

    fout_ = sam_open(ss.str().c_str(), modes[FLAGS_output_flag]); 
    if (!fout_) {
      throw std::runtime_error("Cannot open sam output file");
    }
    sam_hdr_write(fout_, aux->h);
  }
  // start writing to file
  for (int i = 0; i < n_elements; ++i){
    sam_write1(fout_, aux->h, buf[i]); 
    bam_destroy1(buf[i]);
  }
  if (use_file) {
    sam_close(fout_);
    file_id_++;
  }

  VLOG(1) << "Written " << n_elements 
          << " records in "
          << getUs() - start_ts << " us";
}

void SamsPrint::compute() {

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
  const char *modes[] = {"wb", "wbu", "w"};
#else
  FILE* fout;
#endif

  if (use_file) {
    std::stringstream ss;
    ss << out_dir << "/part-"
       << std::setw(6) << std::setfill('0') << file_id;
#ifdef USE_HTSLIB
    fout = sam_open(ss.str().c_str(), modes[FLAGS_output_flag]); 
    sam_hdr_write(fout, aux->h);
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
    // TODO: temporary
    fout_ = sam_open("-", modes[FLAGS_output_flag]); 
    int status = sam_hdr_write(fout, aux->h);
    if (status) {
      LOG(ERROR) << "sam_hdr_write error: " << status;
    }
#else
    fout = stdout;
#endif
  }

  // Buffer to sort output bam
  int      max_bam_records = FLAGS_max_num_records;
  bam1_t** bam_buffer;
  int      bam_buffer_idx = 0;
  if (FLAGS_sort) {
    bam_buffer = (bam1_t**)malloc(FLAGS_max_num_records*sizeof(bam1_t*));
  }

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

        VLOG(1) << "Written batch " << record.start_idx << " to file in "
          << getUs() - start_ts << " us";
      }
    }
    else if (FLAGS_sort) {
      uint64_t start_ts = getUs();
      int batch_num = input.batch_num;
      bseq1_t* seqs = input.seqs;

      for (int i = 0; i < batch_num; i++) {
        if (seqs[i].bams) {   
          for (int j = 0; j < seqs[i].bams->l; j++) {
            bam_buffer[bam_buffer_idx++] = seqs[i].bams->bams[j];
            if (bam_buffer_idx >= max_bam_records) {
              sortAndWriteBamBatch(bam_buffer, max_bam_records, out_dir);
              bam_buffer_idx = 0;
            }
          }  
          free(seqs[i].bams->bams);
          free(seqs[i].bams); seqs[i].bams = NULL;    
        }
      }
      free(seqs);
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
          ss << out_dir << "/part-"
            << std::setw(6) << std::setfill('0') << file_id;

#ifdef USE_HTSLIB
          sam_close(fout);
          fout = sam_open(ss.str().c_str(), modes[FLAGS_output_flag]); 
          sam_hdr_write(fout, aux->h);
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

      VLOG(1) << "Written batch " << input.start_idx << " to file in "
        << getUs() - start_ts << " us";
    }
  }
  if (bam_buffer_idx > 0) {
    sortAndWriteBamBatch(bam_buffer, bam_buffer_idx, out_dir);
  }
#ifdef USE_HTSLIB    
  if (!FLAGS_sort) {
    sam_close(fout);
  }
  bam_hdr_destroy(aux->h);
#else
  if (use_file) {
    fclose(fout);
  }
#endif
}
