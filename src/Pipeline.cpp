#include <boost/smart_ptr.hpp>
#include <boost/thread/thread.hpp>
#include <list>
#include <queue>

#include "bwa/utils.h"
#include "config.h"
#include "Extension.h"
#include "Pipeline.h"  
#include "util.h"

void SeqsProducer::compute() {

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
  std::list<SWRead*>* read_batch = new std::list<SWRead*>;
#endif

  for (int i = 0; i < batch_num; i++) {
    chains[i] = seq2chain(aux, &seqs[i]);
    kv_init(alnreg[i]);

#ifdef OFFLOAD
    SWRead *read_ptr = new SWRead(start_idx, i, aux, 
        seqs+i, chains+i, alnreg+i);

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
  std::list<SWRead*>* new_reads = record.read_batch;

  input_buf[start_idx] = record;

  // push record to output table
  RegionsRecord output;
  output.start_idx = start_idx;
  output.batch_num = batch_num;
  output.seqs = seqs;
  output.alnreg = alnreg;

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
    std::queue<boost::shared_ptr<boost::thread> > stage_workers;

    int task_num[stage_num] = {0};

    // Batch of SWTasks
    ExtParam** task_batch[stage_num];

    for (int i = 0; i < stage_num; i++) {
      task_num[i] = 0;
      task_batch[i] = new ExtParam*[chunk_size];
    }

    // Batch of SWReads
    std::list<SWRead*> read_batch;

    // Table to keep track of each record
    std::unordered_map<uint64_t, int> tasks_remain;
    std::unordered_map<uint64_t, ChainsRecord> input_buf;
    std::unordered_map<uint64_t, RegionsRecord> output_buf;

    // For statistics
    uint64_t swFPGA_time = 0;
    uint64_t extCPU_time = 0;
    int      swFPGA_num  = 0;
    int      extCPU_num  = 0;

    uint64_t wait_time = 0;
    uint64_t nextTask_time = 0;
    uint64_t nextTask_num = 0;

    uint64_t last_batch_ts;
    uint64_t last_output_ts;
    uint64_t last_read_ts;

    uint64_t batch_time  = 0;
    uint64_t output_time = 0;
    int      batch_num   = 0;
    int      read_num    = 0;

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
          uint64_t start_idx;

          ExtParam* param_task;
          start_ts = getUs();
          SWRead::TaskStatus status = (*iter)->nextTask(param_task);
          nextTask_time += getUs() - start_ts; 
          nextTask_num ++;

          int curr_size = 0;

          switch (status) {

            case SWRead::TaskStatus::Successful:

              task_batch[stage_cnt][task_num[stage_cnt]] = param_task;
              task_num[stage_cnt]++;
              if (task_num[stage_cnt] >= chunk_size) {
                VLOG(3) << "nextTask takes " << nextTask_time << " us in " 
                        << nextTask_num << " calls";

                nextTask_time = 0;
                nextTask_num = 0;

                start_ts = getUs();
#ifdef USE_FPGA
                packData(stage_cnt,
                    task_batch[stage_cnt],
                    task_num[stage_cnt],
                    aux->opt);
#endif
                if (!stage_workers.empty()) {
                  stage_workers.front()->join();
                  stage_workers.pop();
                }
                VLOG(3) << "Batch takes " << getUs() - last_batch_ts << " us";
                if (batch_num < 500) {
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
#ifdef USE_FPGA
                boost::shared_ptr<boost::thread> worker(new 
                    boost::thread(&SwFPGA,
                      stage_cnt,
                      task_batch[stage_cnt],
                      task_num[stage_cnt],
                      aux->opt));
                //SwFPGA(task_batch, task_num, aux->opt);
#else
                start_ts = getUs();
                boost::shared_ptr<boost::thread> worker(new 
                    boost::thread(&extendOnCPU,
                      task_batch[stage_cnt],
                      task_num[stage_cnt],
                      aux->opt));
                //extendOnCPU(task_batch, task_num, aux->opt);
#endif
                last_batch_ts = getUs();
                stage_workers.push(worker);

                task_num[stage_cnt] = 0;
                stage_cnt = (stage_cnt + 1) % stage_num;

                swFPGA_time += getUs() - start_ts;
                swFPGA_num ++;
              }
              iter ++;
              break;

            case SWRead::TaskStatus::Pending:
              if (flag_more_reads) {
                // Try to get a new batch
                curr_size = read_batch.size(); 

                start_ts = getUs();
                if (addBatch(read_batch, tasks_remain, input_buf, output_buf)) {
                  iter = read_batch.begin();
                  std::advance(iter, curr_size);
                }
                else {
                  flag_more_reads = false;
                }
                wait_time += getUs() - start_ts;
              }
              else {
                // No more new tasks, must do extend before proceeding
                if (!stage_workers.empty()) {
                  stage_workers.front()->join();
                  stage_workers.pop();
                }
                else {
                  start_ts = getUs();

                  extendOnCPU(task_batch[stage_cnt], task_num[stage_cnt], aux->opt);

                  extCPU_time += getUs() - start_ts;
                  extCPU_num ++;

                  task_num[stage_cnt] = 0;
                  iter++;
                }
              }
              break;

            case SWRead::TaskStatus::Finished:
              // Read is finished, remove from batch
              start_idx = (*iter)->start_idx();
              tasks_remain[start_idx]--;

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

                // Free data in the input record
                freeChains(input_buf[start_idx].chains, 
                    input_buf[start_idx].batch_num);

                tasks_remain.erase(start_idx);
                input_buf.erase(start_idx);
                output_buf.erase(start_idx);

                VLOG(1) << "Produced a region batch in "
                  << getUs() - last_output_ts << " us";
                VLOG(1) << "Currently there are " << tasks_remain.size()
                  << " active batches";

                last_output_ts = getUs();
              }
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

      mem_alnreg_v* alnreg = (mem_alnreg_v*)malloc(batch_num*sizeof(mem_alnreg_v));

      for (int i = 0; i < batch_num; i++) {
        kv_init(alnreg[i]);
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

      freeChains(chains, batch_num);

      RegionsRecord output;
      output.start_idx = record.start_idx;
      output.batch_num = batch_num;
      output.seqs = seqs;
      output.alnreg = alnreg;

      VLOG(1) << "Produced a region batch in "
        << getUs() - last_output_ts << " us";

      pushOutput(output);
    }
  }
}

SeqsRecord RegionsToSam::compute(RegionsRecord const & record) {

  uint64_t start_ts = getUs();
  if (!aux) {
    boost::any var = this->getConst("aux");
    aux = boost::any_cast<ktp_aux_t*>(var);
  }

  int start_idx = record.start_idx;
  int batch_num = record.batch_num;
  mem_alnreg_v* alnreg = record.alnreg;
  bseq1_t* seqs = record.seqs;

  // Post-process each chain before output
  for (int i = 0; i < batch_num; i++) {
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

  SeqsRecord output;
  output.start_idx = start_idx;
  output.batch_num = batch_num;
  output.seqs = seqs;

  VLOG(1) << "Produced a sam batch in "
    << getUs() - start_ts << " us";

  return output;
}

void PrintSam::compute() {
  
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
      }
      freeSeqs(seqs, batch_num);

      // Remove the record from buffer
      record_buf.erase(n_processed);

      n_processed += batch_num;

      VLOG(1) << "Written " << batch_num << " seqs to file in "
        << getUs() - start_ts << " us";
    }
  }
}

