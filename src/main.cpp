#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <glog/logging.h>
#include <string>

#include "bwa/bntseq.h"
#include "bwa/bwa.h"
#include "bwa/bwamem.h"
#include "bwa/kseq.h"
#include "bwa/kvec.h"
#include "bwa/utils.h"

#include "bwa_wrapper.h"
#include "util.h"
#include "FPGAAgent.h"

FPGAAgent* agent;

// global parameters
gzFile fp_idx, fp2_read2 = 0;
void *ko_read1 = 0, *ko_read2 = 0;

void mem_chain2aln_hw(
    ktp_aux_t *aux,
    const bseq1_t *seqs,
    const mem_chain_v* chains,
    mem_alnreg_v *av,
    int batch_num);

int main(int argc, char *argv[]) {
  extern char *bwa_pg;
  extern gzFile fp_idx, fp2_read2;
  extern void *ko_read1, *ko_read2;

  ktp_aux_t aux;
  memset(&aux, 0, sizeof(ktp_aux_t));

  // get the index and the options
  uint64_t start_ts = getUs();
  pre_process(argc-1, argv+1, &aux);
  fprintf(stderr, "Preprocessing time is %dus\n", getUs()-start_ts);

  int batch_num = 0;
  bseq1_t *seqs = bseq_read(1500000, &batch_num, aux.ks, aux.ks2);

  agent = new FPGAAgent("/curr/diwu/prog/acc_lib/bwa-sm/sm-80pe.xclbin", 2000);

  mem_alnreg_v* alnreg = new mem_alnreg_v[batch_num];
  mem_alnreg_v* alnreg_hw = new mem_alnreg_v[batch_num];

  mem_chain_v* chains = new mem_chain_v[batch_num];

  uint64_t cost_sw = 0;
  for (int i = 0; i < batch_num; i++) {
    chains[i] = seq2chain(&aux, &seqs[i]);

    uint64_t start_ts = getUs();
    kv_init(alnreg[i]);

    for (int j = 0; j < chains[i].n; j++) {
      mem_chain_t *p = &chains[i].a[j];

      // call mem_chain2aln to compute baseline
      mem_chain2aln(aux.opt, aux.idx->bns, aux.idx->pac,
          seqs[i].l_seq,
          (uint8_t*)seqs[i].seq,
          p, alnreg+i);
    }
    cost_sw += getUs() - start_ts;
  }
  fprintf(stderr, "Software compute time for %d reads is %dus\n", batch_num, cost_sw);

  start_ts = getUs();
  mem_chain2aln_hw(&aux, seqs, chains, alnreg_hw, batch_num);
  uint64_t cost_hw = getUs() - start_ts;

  fprintf(stderr, "FPGA compute time for %d reads is %dus\n", batch_num, cost_hw);

  regionsCompare(alnreg, alnreg_hw, batch_num);

  // Free the chains
  for (int i = 0; i < batch_num; i++) {
    for (int j = 0; j < chains[i].n; j++) {
      free(chains[i].a[j].seeds);
    }
    free(chains[i].a);
  }
  delete [] chains;

  // Free aligned regions
  for (int i = 0; i < batch_num; i++) {
    free(alnreg[i].a);
    free(alnreg_hw[i].a);
    free(seqs[i].name); 
    free(seqs[i].comment);
    free(seqs[i].seq); 
    free(seqs[i].qual); 
    free(seqs[i].sam);
  }
  free(seqs);
  delete [] alnreg;
  delete [] alnreg_hw;

  delete agent;

  free(aux.opt);
  bwa_idx_destroy(aux.idx);
  kseq_destroy(aux.ks);
  err_gzclose(fp_idx); kclose(ko_read1);

  if (aux.ks2) {
    kseq_destroy(aux.ks2);
    err_gzclose(fp2_read2); kclose(ko_read2);
  }

  return 0;
}
