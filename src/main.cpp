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

// global parameters
gzFile fp_idx, fp2_read2 = 0;
void *ko_read1 = 0, *ko_read2 = 0;


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
  bseq1_t *seqs = bseq_read(150000, &batch_num, aux.ks, aux.ks2);

  smem_aux_t** smems = new smem_aux_t*[batch_num];
  smem_aux_t** smems_new = new smem_aux_t*[batch_num];
  uint64_t cost_sw = 0;
  // convert to 2-bit encoding 
  for (int i = 0; i < batch_num ; i++){
     for (int j =0 ; j < seqs[i].l_seq ;j++){
        seqs[i].seq[j] = seqs[i].seq[j] < 4? seqs[i].seq[j]: nst_nt4_table[(int)seqs[i].seq[j]];
     }
  }
  // the original mem_collect_intv
  for (int i = 0; i < batch_num ; i++){
    smems[i] = smem_aux_init();
    mem_collect_intv (aux.opt, aux.idx->bwt, seqs[i].l_seq,(uint8_t*)seqs[i].seq,smems[i]);
  }
  // the revised mem_collect_intv
  for (int i = 0; i < batch_num; i ++){
    smems_new[i] = smem_aux_init();
    mem_collect_intv_new (aux.opt, aux.idx->bwt, seqs[i].l_seq,(uint8_t*)seqs[i].seq,smems_new[i]);
  }
 
  smemCompare(smems,smems_new,batch_num);
  //regionsCompare(alnreg, alnreg_hw, batch_num);

  // Free aligned regions
  for (int i = 0; i < batch_num; i++) {
    free(seqs[i].name); 
    free(seqs[i].comment);
    free(seqs[i].seq); 
    free(seqs[i].qual); 
    free(seqs[i].sam);
  }
  free(seqs);
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
