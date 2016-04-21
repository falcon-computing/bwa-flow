#include <unordered_map>
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
  int batch_num = seqs_record.batch_num;
  
  mem_chain_v*  chains = (mem_chain_v*)malloc(batch_num*sizeof(mem_chain_v));
  mem_alnreg_v* alnreg = (mem_alnreg_v*)malloc(batch_num*sizeof(mem_alnreg_v));
  std::list<SWRead*>* read_batch = new std::list<SWRead*>;

  for (int i = 0; i < batch_num; i++) {
    chains[i] = seq2chain(aux, &seqs[i]);
    kv_init(alnreg[i]);

    SWRead *read_ptr = new SWRead(i, aux, 
        seqs+i, chains+i, alnreg+i);

    read_batch->push_back(read_ptr); 
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

void ChainsToRegions::compute() {

  boost::any var = this->getConst("aux");
  ktp_aux_t* aux = boost::any_cast<ktp_aux_t*>(var);

  var = this->getConst("chunk_size");
  int chunk_size = boost::any_cast<int>(var);

  // Batch of SWTasks
  ExtParam **task_batch = new ExtParam*[chunk_size];

  // Batch of SWReads
  std::list<SWRead*> read_batch;

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

    bseq1_t* seqs        = record.seqs;
    mem_chain_v* chains  = record.chains;
    int batch_num        = record.batch_num;
    mem_alnreg_v* alnreg = record.alnreg;

#ifndef OFFLOAD
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
    std::list<SWRead*>* read_batch = record.read_batch;
    for (std::list<SWRead*>::iterator iter = read_batch->begin(); 
         iter != read_batch->end();
         iter ++) 
    {
      delete *iter;
    }
    delete read_batch;
    freeChains(chains, batch_num);
#else
    int task_num = 0;
    std::list<SWRead*>* new_reads = record.read_batch;
    
    // copy all new reads to current read_batch
    read_batch.insert(read_batch.end(), new_reads->begin(), new_reads->end());

    delete new_reads;

    DLOG(INFO) << "Add " << read_batch.size() << " new reads to process";

    //std::unordered_map<uint64_t, int> tasks_remain;
    //std::unordered_map<uint64_t, RegionsRecord> output_buf;;

    while (!read_batch.empty()) {

      std::list<SWRead*>::iterator iter = read_batch.begin();
      while (iter != read_batch.end()) {

        uint64_t start_ts;
        ExtParam* param_task;
        SWRead::TaskStatus status = (*iter)->nextTask(param_task);

        switch (status) {

          case SWRead::TaskStatus::Successful:
            task_batch[task_num] = param_task;
            task_num++;
            if (task_num >= chunk_size) {
#ifdef USE_FPGA
              SwFPGA(task_batch, task_num, aux->opt);
#else
              extendOnCPU(task_batch, task_num, aux->opt);
#endif
              task_num = 0;
            }
            ++iter;
            break;

          case SWRead::TaskStatus::Pending:
            // No more tasks, must do extend before proceeding
            extendOnCPU(task_batch, task_num, aux->opt);

            task_num = 0;
            break;

          case SWRead::TaskStatus::Finished:
            //uint64_t start_idx = (*iter)->start_idx();
            //task_remain[start_idx]--;
            //if (task_remain[start_idx] == 0) {
            //  pushOutput(output_buf[start_idx]);
            //}

            // Read is finished, remove from batch
            delete *iter;
            iter = read_batch.erase(iter);
            break;

          default: ;
        }
      }
    }
#endif

    RegionsRecord output;
    output.start_idx = record.start_idx;
    output.batch_num = batch_num;
    output.seqs = seqs;
    output.alnreg = alnreg;

    pushOutput(output);
  }
  delete [] task_batch;
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
