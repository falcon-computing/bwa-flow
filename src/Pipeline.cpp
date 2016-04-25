#include <boost/smart_ptr.hpp>
#include <boost/thread/thread.hpp>
#include <list>
#include <queue>

#include "bwa/utils.h"
#include "Extension.h"
#include "Pipeline.h"  

#define OFFLOAD
#define USE_FPGA

void SeqsProducer::compute() {

  boost::any var = this->getConst("aux");
  ktp_aux_t* aux = boost::any_cast<ktp_aux_t*>(var);

  int num_seqs_produced = 0;
  while (true) {

    // Read from file input, get mem_chains
    int batch_num = 0;
    bseq1_t *seqs = bseq_read(aux->actual_chunk_size, 
        &batch_num, aux->ks, aux->ks2);

    if (!seqs) break;

    LOG(INFO) << "Read " << batch_num << " seqs...";

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

  while (!this->isFinal() && !ready) {
    boost::this_thread::sleep_for(boost::chrono::microseconds(100));
    ready = this->getInput(record);
  }
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
  read_batch.insert(read_batch.end(), new_reads->begin(), new_reads->end());
  DLOG(INFO) << "Add " << new_reads->size() << " new reads to process";

  delete new_reads;

  return true;
}

void ChainsToRegions::compute() {

  boost::any var = this->getConst("aux");
  ktp_aux_t* aux = boost::any_cast<ktp_aux_t*>(var);

#ifdef OFFLOAD
  var = this->getConst("chunk_size");
  int chunk_size = boost::any_cast<int>(var);

  const int stage_num = 2;
  int stage_cnt = 0;
  std::queue<boost::shared_ptr<boost::thread> > stage_workers;

  int task_num[stage_num] = {0};
  blaze::Task_ptr fpga_task[stage_num];

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

  bool flag_need_reads = false;
  bool flag_more_reads = true;

  uint64_t last_batch_ts = blaze::getUs();
  while (flag_more_reads || !read_batch.empty()) { 

    if (read_batch.empty()) {
      // get initial input batch
      if (!addBatch(read_batch, tasks_remain, input_buf, output_buf)) {
        flag_more_reads = false;
      }
    }
    else {
      std::list<SWRead*>::iterator iter = read_batch.begin();

      while (iter != read_batch.end()) {

        uint64_t start_ts;
        uint64_t start_idx;

        ExtParam* param_task;
        SWRead::TaskStatus status = (*iter)->nextTask(param_task);

        int curr_size = 0;

        switch (status) {

          case SWRead::TaskStatus::Successful:

            task_batch[stage_cnt][task_num[stage_cnt]] = param_task;
            task_num[stage_cnt]++;
            if (task_num[stage_cnt] >= chunk_size) {
              start_ts = blaze::getUs();
#ifdef USE_FPGA
              uint64_t pd_ts = blaze::getUs();
              fpga_task[stage_cnt] = packData(task_batch[stage_cnt],
                  task_num[stage_cnt],
                  aux->opt);

              if (!stage_workers.empty()) {
                stage_workers.front()->join();
                stage_workers.pop();
                DLOG(INFO) << "Batch takes " << blaze::getUs() - last_batch_ts << " us";
                DLOG(INFO) << "packData takes " << blaze::getUs() - pd_ts << " us";
              }
              boost::shared_ptr<boost::thread> worker(new 
                  boost::thread(&SwFPGA,
                    task_batch[stage_cnt],
                    fpga_task[stage_cnt], 
                    task_num[stage_cnt],
                    aux->opt));
              //SwFPGA(task_batch, task_num, aux->opt);
#else
              if (!stage_workers.empty()) {
                stage_workers.front()->join();
                stage_workers.pop();
                DLOG(INFO) << "Batch takes " << blaze::getUs() - last_batch_ts << " us";
              }

              boost::shared_ptr<boost::thread> worker(new 
                  boost::thread(&extendOnCPU,
                    task_batch[stage_cnt],
                    task_num[stage_cnt],
                    aux->opt));
              //extendOnCPU(task_batch, task_num, aux->opt);
#endif
              last_batch_ts = blaze::getUs();
              stage_workers.push(worker);

              task_num[stage_cnt] = 0;
              stage_cnt = (stage_cnt + 1) % stage_num;

              swFPGA_time += blaze::getUs() - start_ts;
              swFPGA_num ++;
            }
            iter ++;
            break;

          case SWRead::TaskStatus::Pending:
            if (flag_more_reads) {
              // Try to get a new batch
              curr_size = read_batch.size(); 

              start_ts = blaze::getUs();
              if (addBatch(read_batch, tasks_remain, input_buf, output_buf)) {
                iter = read_batch.begin();
                std::advance(iter, curr_size);
              }
              else {
                flag_more_reads = false;
              }
              wait_time += blaze::getUs() - start_ts;
            }
            else {
              // No more new tasks, must do extend before proceeding
              if (!stage_workers.empty()) {
                stage_workers.front()->join();
                stage_workers.pop();
              }
              else {
                start_ts = blaze::getUs();

                extendOnCPU(task_batch[stage_cnt], task_num[stage_cnt], aux->opt);

                extCPU_time += blaze::getUs() - start_ts;
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

            // Check if corresponding batch is finished
            if (tasks_remain[start_idx] == 0) {
              pushOutput(output_buf[start_idx]);

              // Free data in the input record
              freeChains(input_buf[start_idx].chains, 
                  input_buf[start_idx].batch_num);

              tasks_remain.erase(start_idx);
              input_buf.erase(start_idx);
              output_buf.erase(start_idx);

              DLOG(INFO) << "Pushing output " << start_idx
                         << ", currently there are " << tasks_remain.size()
                         << " active batches.";
            }
            break;

          default: ;
        }
      }
    }
  }
  DLOG(INFO) << swFPGA_num << " batched tasks takes " << swFPGA_time << "us";
  DLOG(INFO) << extCPU_num << " normal tasks takes " << extCPU_time << "us";
  DLOG(INFO) << "Waiting for input takes " << wait_time << "us";

  for (int i = 0; i < stage_num; i++) {
    delete [] task_batch[i];
  }
#else
  while (true) { 
    ChainsRecord record;
    bool ready = this->getInput(record);

    while (!this->isFinal() && !ready) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(100));
      ready = this->getInput(record);
    }
    if (!ready) { 
      // this means isFinal() is true and input queue is empty
      break; 
    }

    bseq1_t* seqs       = record.seqs;
    mem_chain_v* chains = record.chains;
    int batch_num       = record.batch_num;

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
    freeChains(chains, batch_num);

    RegionsRecord output;
    output.start_idx = record.start_idx;
    output.batch_num = batch_num;
    output.seqs = seqs;
    output.alnreg = alnreg;

    pushOutput(output);
  }
#endif
}

SeqsRecord RegionsToSam::compute(RegionsRecord const & record) {

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

    start_ts = kestrelFlow::getUs();
    // Find the next batch in the buffer
    while (record_buf.count(n_processed)) {
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
    }
  }
}

