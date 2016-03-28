#include <stdio.h>
#include <string.h>

#include "kstring.h"
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
#include "structuresnew.h"








// global parameters
ktp_aux_t aux;
gzFile fp_idx, fp2_read2 = 0;
void *ko_read1 = 0, *ko_read2 = 0;
bseq1_t *seqs;
int batch_num = 0;
int64_t n_processed = 0;

int main(int argc, char *argv[])
{
	extern char *bwa_pg;
    extern ktp_aux_t aux;
    extern gzFile fp_idx, fp2_read2 ;
    extern void *ko_read1 , *ko_read2 ;
    extern int batch_num;
	extern int64_t n_processed;

	int i, ret=0;
	double t_real;
	kstring_t pg = {0,0,0};
	t_real = realtime();
	ksprintf(&pg, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
	for (i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]);
    bwa_pg = pg.s;

	if (argc < 2) return usage();
	memset(&aux, 0, sizeof(ktp_aux_t));
	pre_process(argc-1,argv+1);  // get the index and the options


	//seqs = (bseq1_t *)calloc(sizeof(bseq1_t));
	seqs = (bseq1_t *)malloc(sizeof(bseq1_t));
    // initial the datastructures
    smem_aux_t *SMEM;
    mem_chain_v * chain;
    chain = (mem_chain_v *)malloc(sizeof(mem_chain_v));
    mem_alnreg_v *alnreg ;


    while(true)
    {
        batch_num = load_reads();   	// get the reads and stores in seqs ----step0
        if(seqs==0)
            break;
        alnreg = (mem_alnreg_v *) malloc(batch_num*sizeof(mem_alnreg_v));
        SMEM = (smem_aux_t *)malloc(sizeof(smem_aux_t));
        SMEM = smem_aux_init();
        for( i =0; i< batch_num; i++)
        {
            //  from seq to interval: SMEM GENERATION and stores in intervals -----------step1
            seq2intv(&seqs[i],SMEM);                                  //Only one seq at this moment
            // from interval to chain : Suffix array lookup and other operation -------- step2
            intv2chain(SMEM,chain);
            // from chain to regs :Smith Waterman extension --------------step 3
            chain2reg (chain,alnreg);

        }
        // then worker2 and print
    smem_aux_destroy(SMEM);
	//free(SMEM);

    reg2sam (alnreg);
    n_processed += batch_num;

    }

 //   free(hdr_line);
	free(aux.opt);
	bwa_idx_destroy(aux.idx);
	kseq_destroy(aux.ks);
	err_gzclose(fp_idx); kclose(ko_read1);
	if (aux.ks2) {
		kseq_destroy(aux.ks2);
		err_gzclose(fp2_read2); kclose(ko_read2);
	}
	err_fflush(stdout);
	err_fclose(stdout);
	if (ret == 0) {
		fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
		fprintf(stderr, "[%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
	}
	free(bwa_pg);

	return ret;
}
