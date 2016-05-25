#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif
//#include "kstring.h"
#include "bwa/bwamem.h"
#include "bwa/bntseq.h"
#include "bwa/ksw.h"
#include "bwa/kvec.h"
#include "bwa/ksort.h"
#include "bwa/utils.h"
#include "bwa_wrapper.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

mem_chain_v seq2chain(
    ktp_aux_t *aux,
    bseq1_t *seqs
) {
  int i;
  mem_chain_v chn;
  for (i = 0; i < seqs->l_seq; ++i) // convert to 2-bit encoding if we have not done so
    seqs->seq[i] = seqs->seq[i] < 4? seqs->seq[i] : nst_nt4_table[(int)seqs->seq[i]];

  chn = mem_chain(aux->opt, aux->idx->bwt, aux->idx->bns, seqs->l_seq, (uint8_t*)seqs->seq, 0);         // the 0 should be reconsidered
  chn.n = mem_chain_flt(aux->opt, chn.n, chn.a);
  mem_flt_chained_seeds(aux->opt, aux->idx->bns, aux->idx->pac, seqs->l_seq, (uint8_t*)seqs->seq, chn.n, chn.a);

  return chn;
}
/*
void reg2sam(
    ktp_aux_t* aux,
    bseq1_t* seqs,
    int batch_num,
    int64_t n_processed,
    mem_alnreg_v* alnreg
) {
  // the operation between worker1 and worker2
  mem_pestat_t pes[4];
  mem_pestat(aux->opt, aux->idx->bns->l_pac, batch_num, alnreg, pes);
  for (int i = 0; i < batch_num/2; i++) {
    mem_sam_pe(
        aux->opt,
        aux->idx->bns,
        aux->idx->pac,
        pes,
        (n_processed>>1)+i,
        &seqs[i<<1],
        &alnreg[i<<1]);
  }
  // print
  for (int i = 0; i < batch_num; ++i) {
    if (seqs[i].sam) err_fputs(seqs[i].sam, stdout);
  }
}
*/
void freeChains(mem_chain_v* chains, int batch_num) {
  // Free the chains
  for (int i = 0; i < batch_num; i++) {
    for (int j = 0; j < chains[i].n; j++) {
      free(chains[i].a[j].seeds);
    }
    free(chains[i].a);
  }
  free(chains);
}

void freeAligns(mem_alnreg_v* alnreg, int batch_num) {
  for (int i = 0; i < batch_num; i++) {
    free(alnreg[i].a); 
  }
  free(alnreg);
}

/*
void freeSeqs(bseq1_t* seqs, int batch_num) {
  for (int i = 0; i < batch_num; i++) {
    free(seqs[i].name); 
    free(seqs[i].comment);
    free(seqs[i].seq); 
    free(seqs[i].qual); 
    free(seqs[i].sam);
  }
  free(seqs);
}
*/
