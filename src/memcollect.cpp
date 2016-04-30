#include <stdlib.h>
#include "bwa_wrapper.h"
#include "bwa/bwamem.h"

smem_aux_t *smem_aux_init()
{
	smem_aux_t *a;
	a = (smem_aux_t *)calloc(1, sizeof(smem_aux_t));
	a->tmpv[0] = (bwtintv_v*)calloc(1, sizeof(bwtintv_v));
	a->tmpv[1] = (bwtintv_v*)calloc(1, sizeof(bwtintv_v));
	return a;
}


 void mem_collect_intv_new(const mem_opt_t *opt, const bwt_t *bwt, int len, const uint8_t *seq, smem_aux_t *a)
{

   mem_collect_intv(opt,bwt,len,seq,a);

}
