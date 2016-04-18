#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include "bwa/bntseq.h"
#include "bwa/bwa.h"
#include "bwa/bwamem.h"
#include "bwa/kseq.h"
#include "bwa/kvec.h"
#include "bwa/utils.h"

#include "bwa_wrapper.h"
//#include "util.h"
//#include <glog/logging.h>
#include <string>

// global parameters
gzFile fp_idx, fp2_read2 = 0;
void *ko_read1 = 0, *ko_read2 = 0;

/*
void mem_chain2aln_hw(
    ktp_aux_t *aux,
    const bseq1_t *seqs,
    const MemChainVector* chains,
    mem_alnreg_v *av,
    int batch_num);
*/

int main(int argc, char *argv[]) {

	double t_real = realtime();

  // Preprocessing
  // TODO: move these to the global var lists
  extern char *bwa_pg;
  extern gzFile fp_idx, fp2_read2;
  extern void *ko_read1, *ko_read2;
  ktp_aux_t aux;

  memset(&aux, 0, sizeof(ktp_aux_t));

  kstring_t pg = {0,0,0};
  ksprintf(&pg, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
  for (int i = 1; i < argc; i++) {
    ksprintf(&pg, " %s", argv[i]);
  }
  bwa_pg = pg.s;

  // get the index and the options
  pre_process(argc-1, argv+1, &aux);

  // Read from file input, get mem_chains
  uint64_t n_processed = 0;
  while (true) {
	  double t_real = realtime();

    int batch_num = 0;
    bseq1_t *seqs = bseq_read(aux.actual_chunk_size, 
        &batch_num, aux.ks, aux.ks2);

		fprintf(stderr, "[M::%s] read %d sequences...\n", __func__, batch_num);

    if (!seqs) break;

    mem_chain_v* chains = (mem_chain_v*)malloc(batch_num*sizeof(mem_chain_v));

    for (int i = 0; i < batch_num; i++) {
      chains[i] = seq2chain(&aux, &seqs[i]);
    }

    // Process each chain to get the aligned regions
    mem_alnreg_v* alnreg = (mem_alnreg_v*)malloc(batch_num*sizeof(mem_alnreg_v));

    for (int i = 0; i < batch_num; i++) {
      kv_init(alnreg[i]);
      for (int j = 0; j < chains[i].n; j++) {
        mem_chain2aln(
            aux.opt, 
            aux.idx->bns, 
            aux.idx->pac,
            seqs[i].l_seq,
            (uint8_t*)seqs[i].seq,
            &chains[i].a[j],
            alnreg+i);
      }
    }
    freeChains(chains, batch_num);

    // Post-process each chain and output to sam
    for (int i = 0; i < batch_num; i++) {
      alnreg[i].n = mem_sort_dedup_patch(
          aux.opt, 
          aux.idx->bns, 
          aux.idx->pac, 
          (uint8_t*)seqs[i].seq, 
          alnreg[i].n, 
          alnreg[i].a);

      for (int j = 0; j < alnreg[i].n; j++) {
        mem_alnreg_t *p = &alnreg[i].a[j];
        if (p->rid >= 0 && aux.idx->bns->anns[p->rid].is_alt)
          p->is_alt = 1;
      }
    }
    reg2sam(&aux, seqs, batch_num, n_processed, alnreg);
    freeAligns(alnreg, batch_num);
    freeSeqs(seqs, batch_num);

    n_processed += batch_num;

    fprintf(stderr, "[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
  }
  
  // Free all global variables
  //delete agent;
  free(aux.opt);
  bwa_idx_destroy(aux.idx);
  kseq_destroy(aux.ks);
  err_gzclose(fp_idx); 
  kclose(ko_read1);

  if (aux.ks2) {
    kseq_destroy(aux.ks2);
    err_gzclose(fp2_read2); kclose(ko_read2);
  }

	err_fflush(stdout);
	err_fclose(stdout);

  fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
  fprintf(stderr, "[%s] CMD:", __func__);
  for (int i = 0; i < argc; ++i) {
    fprintf(stderr, " %s", argv[i]);
  }
  fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());

	free(bwa_pg);
  return 0;
}
