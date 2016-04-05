#include <stdio.h>
#include <string.h>

//#include "kstring.h"
#include "utils.h"


#include <zlib.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "bntseq.h"
#include "kseq.h"
#include "bwa_wrapper.h"

// global parameters

gzFile fp_idx, fp2_read2 = 0;
void *ko_read1 = 0, *ko_read2 = 0;

void mem_chain2aln_hw(
    ktp_aux_t *aux,
    const bseq1_t *seqs, 
    const MemChainVector* chains, 
    mem_alnreg_v *av,
    int batch_num);

int main(int argc, char *argv[])
{
  extern char *bwa_pg;
  extern gzFile fp_idx, fp2_read2 ;
  extern void *ko_read1 , *ko_read2 ;

  ktp_aux_t aux;
  memset(&aux, 0, sizeof(ktp_aux_t));

  pre_process(argc-1,argv+1,&aux);  // get the index and the options

  mem_alnreg_v *alnreg ;
  mem_alnreg_v *alnreg_hw ;

  int batch_num = 0;
  bseq1_t *seqs = bseq_read(150000,&batch_num,aux.ks,aux.ks2);

  alnreg = (mem_alnreg_v *) malloc(batch_num*sizeof(mem_alnreg_v));
  alnreg_hw = (mem_alnreg_v *) malloc(batch_num*sizeof(mem_alnreg_v));

  MemChainVector* chains = new MemChainVector[batch_num];

  for(int i=0; i<batch_num; i++)
  {
    chains[i] = seq2chain(&aux, &seqs[i]);
    chains[i].id_read = i;

    kv_init(alnreg[i]);
    for (int j=0; j<chains[i].n; j++) {
      mem_chain_t *p = &chains[i].a[j];

      // call mem_chain2aln to compute baseline
		  mem_chain2aln( aux.opt, aux.idx->bns, aux.idx->pac, 
          seqs[i].l_seq, (uint8_t*)seqs[i].seq, p, alnreg+i);
    }
  }
  mem_chain2aln_hw(&aux, seqs, chains, alnreg_hw, batch_num);

  // free the chains
  for(int i=0; i<batch_num; i++) {
    for (int j=0; j<chains[i].n; j++) {
		  free(chains[i].a[j].seeds);
    }
    free(chains[i].a);
  }
  delete [] chains;
  
  reg_dump(alnreg, alnreg_hw, batch_num);

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
