#include <boost/asio.hpp>
#include <boost/function_types/result_type.hpp>
#include <boost/make_shared.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread.hpp>
#include <list>
#include <queue>

#include "mpi.h"

#include "bwa/utils.h"
#include "config.h"
#include "Extension.h"
#include "Pipeline.h"  
#include "util.h"

// Helper functions in obj serialization/deserialization
template <typename T>
static inline void putT(std::stringstream &ss, T value) {
  ss.write(reinterpret_cast<char*>(&value), sizeof(value));
}

template <typename T>
static inline void getT(std::stringstream &ss, T &value) {
  ss.read(reinterpret_cast<char*>(&value), sizeof(value));
}

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

static inline void getStr(std::stringstream &ss, char* &dst) {
  size_t length = 0;
  getT(ss, length);

  if (length > 0) {
    dst = (char*)malloc(length+1);
    ss.read(dst, length);
    dst[length] = '\0';
  }
}

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

  boost::any var = this->getConst("aux");
  ktp_aux_t* aux = boost::any_cast<ktp_aux_t*>(var);

  var = this->getConst("mpi_rank");
  int rank = boost::any_cast<int>(var);

  var = this->getConst("mpi_nprocs");
  int nprocs = boost::any_cast<int>(var);

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

      // First query a process to send data to
      int proc_id = -1;
      bwaMPIRecv(&proc_id, 1, MPI::INT, MPI::ANY_SOURCE, SEQ_DP_QUERY);

      bwaMPISend(&length, 1, MPI::INT, proc_id, SEQ_DP_LENGTH);

      bwaMPISend(ser_data.c_str(), length, MPI::CHAR, proc_id, SEQ_DP_DATA);

      VLOG(2) << "Sending batch " << record.start_idx
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

  boost::any var = this->getConst("mpi_rank");
  int rank = boost::any_cast<int>(var);

  var = this->getConst("mpi_nprocs");
  int nprocs = boost::any_cast<int>(var);

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

void SeqsRead::compute() {

  boost::any var = this->getConst("aux");
  ktp_aux_t* aux = boost::any_cast<ktp_aux_t*>(var);

  int num_seqs_produced = 0;

  while (true) {
    uint64_t start_ts = getUs();

    // Read from file input, get mem_chains
    int batch_num = 0;
    bseq1_t *seqs = bseq_read(aux->actual_chunk_size, 
        &batch_num, aux->ks, aux->ks2);

    if (!seqs) break;

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

  uint64_t start_ts = getUs();

  if (!aux) {
    boost::any var = this->getConst("aux");
    aux = boost::any_cast<ktp_aux_t*>(var);
  }

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

ChainsRecord SeqsToChains::compute(SeqsRecord const & seqs_record) {

  uint64_t start_ts = getUs();

  if (!aux) {
    boost::any var = this->getConst("aux");
    aux = boost::any_cast<ktp_aux_t*>(var);
  }

  bseq1_t* seqs = seqs_record.seqs;
  int start_idx = seqs_record.start_idx;
  int batch_num = seqs_record.batch_num;

  mem_chain_v*  chains = (mem_chain_v*)malloc(batch_num*sizeof(mem_chain_v));
  mem_alnreg_v* alnreg = (mem_alnreg_v*)malloc(batch_num*sizeof(mem_alnreg_v));

#ifdef OFFLOAD
  std::vector<int>* chains_idxes = new std::vector<int>[batch_num];
  std::list<SWRead*>* read_batch = new std::list<SWRead*>;
#endif

  for (int i = 0; i < batch_num; i++) {
    chains[i] = seq2chain(aux, &seqs[i]);
    kv_init(alnreg[i]);

#ifdef OFFLOAD
    SWRead *read_ptr = new SWRead(start_idx, i, aux, 
        seqs+i, chains+i, alnreg+i, chains_idxes+i);

    read_batch->push_back(read_ptr); 
#endif
  }

  ChainsRecord ret;
  ret.start_idx = seqs_record.start_idx;
  ret.batch_num = batch_num;
  ret.seqs = seqs;
  ret.chains = chains;
  ret.alnreg = alnreg;
#ifdef OFFLOAD
  ret.chains_idxes = chains_idxes;
  ret.read_batch = read_batch;
#endif

  VLOG(1) << "Produced a chain batch in "
    << getUs() - start_ts << " us";

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
  //read_batch.insert(read_batch.end(), new_reads->begin(), new_reads->end());
  read_batch.splice(read_batch.end(), *new_reads);
  VLOG(2) << "Add " << new_reads->size() << " new reads to process";

  delete new_reads;

  return true;
}

void ChainsToRegions::compute(int wid) {

  boost::any var = this->getConst("aux");
  ktp_aux_t* aux = boost::any_cast<ktp_aux_t*>(var);

#ifdef OFFLOAD
  if (wid == 0) {
#else
  if (wid == -1) {
#endif
    VLOG(1) << "Worker " << wid << " is working on FPGA";
    var = this->getConst("chunk_size");
    int chunk_size = boost::any_cast<int>(var);

    const int stage_num = 2;
    int stage_cnt = 0;

    // Create io_service to post work to the group
    boost::asio::io_service ios;
    boost::asio::io_service::work ios_work(ios);

    // Create a thread group
    boost::thread_group stage_workers;
    stage_workers.create_thread(
        boost::bind(&boost::asio::io_service::run, &ios));
    
    boost::unique_future<void> extension_event;

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
          start_ts = getUs();
          SWRead::TaskStatus status = (*iter)->nextTask(param_task);
          nextTask_time += getUs() - start_ts; 
          nextTask_num ++;

          int curr_size = 0;

          switch (status) {

            case SWRead::TaskStatus::Successful:

              task_batch[stage_cnt][task_num] = param_task;
              task_num++;
              if (task_num >= chunk_size) {
                VLOG(3) << "nextTask takes " << nextTask_time << " us in " 
                        << nextTask_num << " calls";

                nextTask_time = 0;
                nextTask_num = 0;

#ifdef USE_FPGA
                packData(stage_cnt,
                    task_batch[stage_cnt],
                    task_num,
                    aux->opt);
#endif
                if (extension_event.valid()) {
                  uint64_t start_ts = getUs();
                  extension_event.wait();

                  VLOG(3) << "Waiting for last batch to finish takes "
                    << getUs() - start_ts << " us";
                }

                VLOG(3) << "Process SWRead::Pending takes " << pending_time << " us";
                VLOG(3) << "Process SWRead::Finish takes " << finish_time << " us";
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

                typedef boost::packaged_task<void> task_t;
#ifdef USE_FPGA
                boost::shared_ptr<task_t> extension_task =
                    boost::make_shared<task_t>(boost::bind(&SwFPGA,
                        stage_cnt,
                        task_batch[stage_cnt],
                        task_num,
                        aux->opt));
                //SwFPGA(task_batch, task_num, aux->opt);
#else
                boost::shared_ptr<task_t> extension_task =
                    boost::make_shared<task_t>(boost::bind(&extendOnCPU,
                        task_batch[stage_cnt],
                        task_num,
                        aux->opt));
                //extendOnCPU(task_batch, task_num, aux->opt);
#endif
                extension_event = extension_task->get_future();
                ios.post(boost::bind(&task_t::operator(), extension_task));

                last_batch_ts = getUs();

                task_num = 0;
                stage_cnt = (stage_cnt + 1) % stage_num;
              }
              iter ++;
              break;

            case SWRead::TaskStatus::Pending:
              start_ts = getUs();
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
                if (extension_event.valid() && !extension_event.is_ready()) {
                  uint64_t start_ts = getUs();
                  extension_event.wait();
                  VLOG(3) << "Waiting for pending tasks to finish takes "
                    << getUs() - start_ts << " us";
                }
                else {
                  extendOnCPU(task_batch[stage_cnt], task_num, aux->opt);

                  task_num = 0;
                  iter++;
                }
              }
              pending_time += getUs() - start_ts;
              break;

            case SWRead::TaskStatus::Finished:
              start_ts = getUs();
              // Read is finished, remove from batch
              start_idx = (*iter)->start_idx();
              tasks_remain[start_idx]--;

              delete *iter;

              iter = read_batch.erase(iter);

              // Collecting throughput for reads
              if (read_num < 1000) {
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

                // Free data in the input record
                //freeChains(input_buf[start_idx].chains, 
                //    input_buf[start_idx].batch_num);

                tasks_remain.erase(start_idx);
                input_buf.erase(start_idx);
                output_buf.erase(start_idx);

                VLOG(1) << "Produced a region batch in "
                  << getUs() - last_output_ts << " us";
                VLOG(1) << "Currently there are " << tasks_remain.size()
                  << " active batches";

                last_output_ts = getUs();
              }
              finish_time += getUs() - start_ts;
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
  else {
    VLOG(1) << "Worker " << wid << " is working on CPU";

    uint64_t last_output_ts;
    uint64_t last_read_ts;

    while (true) { 
      ChainsRecord record;
      bool ready = this->getInput(record);

      uint64_t start_ts = getUs();
      while (!this->isFinal() && !ready) {
        boost::this_thread::sleep_for(boost::chrono::microseconds(100));
        ready = this->getInput(record);
      }
      if (!ready) { 
        // this means isFinal() is true and input queue is empty
        break; 
      }
      VLOG(2) << "Wait for input takes " << getUs() - start_ts << " us";

      last_output_ts = getUs();
      last_read_ts = getUs();

      bseq1_t* seqs       = record.seqs;
      mem_chain_v* chains = record.chains;
      int batch_num       = record.batch_num;

#ifdef OFFLOAD
      // Free all read batches
      std::list<SWRead*>* read_batch = record.read_batch;
      for (std::list<SWRead*>::iterator iter = read_batch->begin();
          iter != read_batch->end();
          iter ++)
      {
        delete *iter;
      }
      delete read_batch;
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

      VLOG(1) << "Produced a region batch in "
        << getUs() - last_output_ts << " us";

      pushOutput(output);
    }
  }
}

SeqsRecord RegionsToSam::compute(RegionsRecord const & record) {

  uint64_t start_ts = getUs();
  uint64_t seedcov_time = 0;

  if (!aux) {
    boost::any var = this->getConst("aux");
    aux = boost::any_cast<ktp_aux_t*>(var);
  }

  int start_idx        = record.start_idx;
  int batch_num        = record.batch_num;
  mem_chain_v* chains  = record.chains;
  mem_alnreg_v* alnreg = record.alnreg;
  bseq1_t* seqs        = record.seqs;
  std::vector<int>* chains_idxes = record.chains_idxes;

  for (int i = 0; i < batch_num; i++) {
#ifdef OFFLOAD
    uint64_t start_ts = getUs();
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
    seedcov_time += getUs() - start_ts;
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
#ifdef OFFLOAD
  if (chains_idxes) {
    delete [] chains_idxes;
    freeChains(chains, batch_num);
  }
#endif
  VLOG(2) << "Seed coverage time is " << seedcov_time << " us";

  mem_pestat_t pes[4];
  mem_pestat(aux->opt, aux->idx->bns->l_pac, batch_num, alnreg, pes);
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

void SamsSend::compute() {

  boost::any var = this->getConst("mpi_rank");
  int rank = boost::any_cast<int>(var);

  var = this->getConst("mpi_nprocs");
  int nprocs = boost::any_cast<int>(var);

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
      // Serialize data and send to master
      
      std::string ser_data = serialize(&input);

      int length = ser_data.length();

      if (length <= 0) {
        throw std::runtime_error("Possible overflow of msg length");
      }
      // Send proc_id to master to let master receive following msg
      bwaMPISend(&rank, 1, MPI::INT, MASTER_RANK, SAM_RV_QUERY);

      uint64_t start_ts = getUs();
      bwaMPISend(&length, 1, MPI::INT, MASTER_RANK, SAM_RV_LENGTH);

      bwaMPISend(ser_data.c_str(), length,
          MPI::CHAR, MASTER_RANK, SAM_RV_DATA);

      VLOG(2) << "Sending batch " << input.start_idx
        << " to master takes " << getUs() - start_ts << " us";

      //freeSeqs(input.seqs, input.batch_num);
      for (int i = 0; i < input.batch_num; i++) {
        free(input.seqs[i].sam);
      }
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
  
  boost::any var = this->getConst("aux");
  ktp_aux_t* aux = boost::any_cast<ktp_aux_t*>(var);

  var = this->getConst("mpi_nprocs");
  int nprocs = boost::any_cast<int>(var);

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
  
  getT(ss, start_idx);
  getT(ss, batch_num);

  bseq1_t* seqs = (bseq1_t*)malloc(batch_num*sizeof(bseq1_t));

  for (int i = 0; i < batch_num; i++) {
    bseq1_t* seq = &seqs[i];
    memset(seq, 0, sizeof(bseq1_t));
    
    getStr(ss, seq->sam);
  }

  SeqsRecord output;
  output.start_idx = start_idx;
  output.batch_num = batch_num;
  output.seqs = seqs;

  return output;
}

void SamsPrint::compute() {
  
  boost::any var = this->getConst("aux");
  ktp_aux_t* aux = boost::any_cast<ktp_aux_t*>(var);

  uint64_t n_processed = 0;

  std::unordered_map<uint64_t, SeqsRecord> record_buf;
  // NOTE: input may be out-of-order
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

    // Add the current input to buffer first
    record_buf[input.start_idx] = input;

    // Find the next batch in the buffer
    while (record_buf.count(n_processed)) {
      start_ts = getUs();
      
      SeqsRecord record = record_buf[n_processed];

      int      batch_num = record.batch_num;
      bseq1_t* seqs      = record.seqs;

      //reg2sam(aux, curr_seqs, curr_batch_num, n_processed, curr_alnreg);
      for (int i = 0; i < batch_num; ++i) {
        if (seqs[i].sam) err_fputs(seqs[i].sam, stdout);
        free(seqs[i].sam);
      }
      //freeSeqs(seqs, batch_num);
      free(seqs);

      // Remove the record from buffer
      record_buf.erase(n_processed);

      n_processed += batch_num;

      VLOG(1) << "Written " << batch_num << " seqs to file in "
        << getUs() - start_ts << " us";
    }
  }
}
