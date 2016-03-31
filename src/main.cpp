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


int main(int argc, char *argv[])
{
	extern char *bwa_pg;
    extern gzFile fp_idx, fp2_read2 ;
    extern void *ko_read1 , *ko_read2 ;
    int batch_num = 0;
    int i = 0;
   	ktp_aux_t aux;
	bseq1_t *seqs;
	memset(&aux, 0, sizeof(ktp_aux_t));
	pre_process(argc-1,argv+1,&aux);  // get the index and the options
	seqs = (bseq1_t *)malloc(sizeof(bseq1_t));
    MemChainVector chain;
    mem_alnreg_v *alnreg ;
    mem_alnreg_v *alnreg_hw ;

        seqs = bseq_read(150000,&batch_num,aux.ks,aux.ks2);
        alnreg = (mem_alnreg_v *) malloc(batch_num*sizeof(mem_alnreg_v));
        alnreg_hw = (mem_alnreg_v *) malloc(batch_num*sizeof(mem_alnreg_v));
        for(i = 0; i < batch_num; i++)
            {
                chain = seq2chain(&aux, &seqs[i]);
                chain2reg_testhw(&aux,seqs,chain,&alnreg[i],&alnreg_hw[i]);
            }

        reg_dump(alnreg,alnreg_hw,batch_num);



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
