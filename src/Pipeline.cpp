#include "bwa/utils.h"

#include "Pipeline.h"  
#include "Extension.h"

#define OFFLOAD

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
  ret.read_batch = read_batch;

  return ret;
}

bool ChainsToRegions::addBatch() {
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

  input_buf_[start_idx] = record;

  // push record to output table
  RegionsRecord output;
  output.start_idx = start_idx;
  output.batch_num = batch_num;
  output.seqs = seqs;
  output.alnreg = alnreg;

  output_buf_[start_idx] = output;
  tasks_remain_[start_idx] = batch_num;

  // copy all new reads to current read_batch
  read_batch_.insert(read_batch_.end(), new_reads->begin(), new_reads->end());
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

  // Batch of SWTasks
  int task_num = 0;

  task_batch_ = new ExtParam*[chunk_size];

  // For statistics
  uint64_t swFPGA_time = 0;
  uint64_t extCPU_time = 0;
  int      swFPGA_num  = 0;
  int      extCPU_num  = 0;

  bool flag_need_reads = false;
  bool flag_more_reads = true;

  while (flag_more_reads || !read_batch_.empty()) { 

    if (read_batch_.empty()) {
      // get initial input batch
      std::list<SWRead*>::iterator iter;
      if (!addBatch()) {
        flag_more_reads = false;
      }
    }
    else {
      std::list<SWRead*>::iterator iter = read_batch_.begin();

      while (iter != read_batch_.end()) {

        uint64_t start_ts;
        uint64_t start_idx;

        ExtParam* param_task;
        SWRead::TaskStatus status = (*iter)->nextTask(param_task);

        int curr_size = 0;

        switch (status) {

          case SWRead::TaskStatus::Successful:
            task_batch_[task_num] = param_task;
            task_num++;
            if (task_num >= chunk_size) {
              start_ts = blaze::getUs();
#ifdef USE_FPGA
              SwFPGA(task_batch_, task_num, aux->opt);
#else
              extendOnCPU(task_batch_, task_num, aux->opt);
#endif
              iter = read_batch_.begin() ;  // go back to the start

              swFPGA_time += blaze::getUs() - start_ts;
              swFPGA_num ++;

              task_num = 0;
            }
            else {
              ++iter;
            }
            break;

          case SWRead::TaskStatus::Pending:
            if (flag_more_reads) {
              // Try to get a new batch
              curr_size = read_batch_.size(); 

              if (addBatch()) {
                iter = read_batch_.begin();
                std::advance(iter, curr_size);
              }
              else {
                flag_more_reads = false;
              }
            }
            else {
              // No more new tasks, must do extend before proceeding
              start_ts = blaze::getUs();

              extendOnCPU(task_batch_, task_num, aux->opt);

              extCPU_time += blaze::getUs() - start_ts;
              extCPU_num ++;

              iter = read_batch_.begin() ;  // go back to the start

              task_num = 0;
            }
            break;

          case SWRead::TaskStatus::Finished:
            // Read is finished, remove from batch
            start_idx = (*iter)->start_idx();

            delete *iter;

            iter = read_batch_.erase(iter);
            tasks_remain_[start_idx]--;

            // Check if corresponding batch is finished
            if (tasks_remain_[start_idx] == 0) {
              pushOutput(output_buf_[start_idx]);

              // Free data in the input record
              freeChains(input_buf_[start_idx].chains, 
                  input_buf_[start_idx].batch_num);

              tasks_remain_.erase(start_idx);
              input_buf_.erase(start_idx);
              output_buf_.erase(start_idx);

              DLOG(INFO) << "Pushing output " << start_idx
                         << ", currently there are " << tasks_remain_.size()
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

  delete [] task_batch_;
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

