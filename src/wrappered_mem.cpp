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
#include "bwamem.h"
#include "bntseq.h"
#include "ksw.h"
#include "kvec.h"
#include "ksort.h"
#include "utils.h"
#include "bwa_wrapper.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

void seq2intv(ktp_aux_t *aux,bseq1_t *seqs,smem_aux_t *SMEM)
{
    int i ;
  	for (i = 0; i < seqs->l_seq; ++i) // convert to 2-bit encoding if we have not done so
		seqs->seq[i] = seqs->seq[i] < 4? seqs->seq[i] : nst_nt4_table[(int)seqs->seq[i]];
    mem_collect_intv(aux->opt,aux->idx->bwt,seqs->l_seq,(uint8_t *)seqs->seq,SMEM);
    SMEM->id_read = seqs->id;
}


MemChainVector seq2chain(ktp_aux_t *aux,bseq1_t *seqs)
{
	int i;
	MemChainVector chain_wid;
	mem_chain_v chn;
	for (i = 0; i < seqs->l_seq; ++i) // convert to 2-bit encoding if we have not done so
		seqs->seq[i] = seqs->seq[i] < 4? seqs->seq[i] : nst_nt4_table[(int)seqs->seq[i]];

	chn = mem_chain(aux->opt, aux->idx->bwt, aux->idx->bns, seqs->l_seq, (uint8_t*)seqs->seq, 0);         // the 0 should be reconsidered
	chn.n = mem_chain_flt(aux->opt, chn.n, chn.a);
	mem_flt_chained_seeds(aux->opt, aux->idx->bns, aux->idx->pac, seqs->l_seq, (uint8_t*)seqs->seq, chn.n, chn.a);
	chain_wid.id_read = seqs->id;
	chain_wid.n = chn.n;
	chain_wid.m = chn.m;
	chain_wid.a = chn.a;
    return chain_wid;
}


void chain2reg(ktp_aux_t *aux,bseq1_t *seqs,MemChainVector chn,mem_alnreg_v *alnreg)
{
    int id = chn.id_read;
    int i;
	kv_init(*alnreg);
	for (i = 0; i < chn.n; ++i) {
		mem_chain_t *p = &chn.a[i];
		if (bwa_verbose >= 4) err_printf("* ---> Processing chain(%d) <---\n", i);
		mem_chain2aln(aux->opt, aux->idx->bns, aux->idx->pac, seqs[id].l_seq, (uint8_t*)seqs[id].seq, p, alnreg);
		free(chn.a[i].seeds);
	}
	free(chn.a);
    alnreg->n = mem_sort_dedup_patch(aux->opt, aux->idx->bns, aux->idx->pac, (uint8_t*)seqs[id].seq, alnreg->n, alnreg->a);
	if (bwa_verbose >= 4) {
		err_printf("* %ld chains remain after removing duplicated chains\n", alnreg->n);
		for (i = 0; i < alnreg->n; ++i) {
			mem_alnreg_t *p = &alnreg->a[i];
			printf("** %d, [%d,%d) <=> [%ld,%ld)\n", p->score, p->qb, p->qe, (long)p->rb, (long)p->re);
		}
	}
	for (i = 0; i < alnreg->n; ++i) {
		mem_alnreg_t *p = &alnreg->a[i];
		if (p->rid >= 0 && aux->idx->bns->anns[p->rid].is_alt)
			p->is_alt = 1;
	}
}


void reg2sam(ktp_aux_t *aux,bseq1_t *seqs,int batch_num,int64_t n_processed,mem_alnreg_v *alnreg)
{
    int i ;
    mem_pestat_t pes[4];
    mem_pestat(aux->opt, aux->idx->bns->l_pac, batch_num, alnreg, pes);               // the operation between worker1 and worker2
    for (i=0; i<batch_num/2; i++)
        {
            mem_sam_pe(aux->opt,aux->idx->bns,aux->idx->pac,pes,n_processed+i,&seqs[i<<1],&alnreg[i<<1]);
            free(alnreg[i<<1|0].a); free(alnreg[i<<1|1].a);
        }
    free(alnreg);
        // print
    for (i = 0; i < batch_num; ++i)
        {
			if (seqs[i].sam) err_fputs(seqs[i].sam, stdout);
			free(seqs[i].name); free(seqs[i].comment);
			free(seqs[i].seq); free(seqs[i].qual); free(seqs[i].sam);
        }
    free(seqs);

}

