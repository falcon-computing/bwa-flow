#include <unordered_map>
#include "bwa/utils.h"

#include "Pipeline.h"  

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
  
  mem_chain_v* chains = (mem_chain_v*)malloc(batch_num*sizeof(mem_chain_v));
  for (int i = 0; i < batch_num; i++) {
    chains[i] = seq2chain(aux, &seqs[i]);
  }

  ChainsRecord ret;
  ret.start_idx = seqs_record.start_idx;
  ret.batch_num = batch_num;
  ret.seqs = seqs;
  ret.chains = chains;

  return ret;
}

void ChainsToRegions::compute() {

  boost::any var = this->getConst("aux");
  ktp_aux_t* aux = boost::any_cast<ktp_aux_t*>(var);

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

