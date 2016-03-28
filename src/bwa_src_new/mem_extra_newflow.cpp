#include <limits.h>
#include "bwa.h"
#include "bwamem.h"
#include "bntseq.h"
#include "kstring.h"

/***************************
 * SMEM iterator interface *
 ***************************/

struct __smem_i {
	const bwt_t *bwt;
	const uint8_t *query;
	int start, len;
	int min_intv, max_len;
	uint64_t max_intv;
	bwtintv_v *matches; // matches; to be returned by smem_next()
	bwtintv_v *sub;     // sub-matches inside the longest match; temporary
	bwtintv_v *tmpvec[2]; // temporary arrays
};

smem_i *smem_itr_init(const bwt_t *bwt)
{
	smem_i *itr;
	itr = (smem_i *)calloc(1, sizeof(smem_i));
	itr->bwt = bwt;
	itr->tmpvec[0] = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));
	itr->tmpvec[1] = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));
	itr->matches   = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));
	itr->sub       = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));
	itr->min_intv = 1;
	itr->max_len  = INT_MAX;
	itr->max_intv = 0;
	return itr;
}

void smem_itr_destroy(smem_i *itr)
{
	free(itr->tmpvec[0]->a); free(itr->tmpvec[0]);
	free(itr->tmpvec[1]->a); free(itr->tmpvec[1]);
	free(itr->matches->a);   free(itr->matches);
	free(itr->sub->a);       free(itr->sub);
	free(itr);
}

void smem_set_query(smem_i *itr, int len, const uint8_t *query)
{
	itr->query = query;
	itr->start = 0;
	itr->len = len;
}

void smem_config(smem_i *itr, int min_intv, int max_len, uint64_t max_intv)
{
	itr->min_intv = min_intv;
	itr->max_len  = max_len;
	itr->max_intv = max_intv;
}

const bwtintv_v *smem_next(smem_i *itr)
{
	int ori_start;
	itr->tmpvec[0]->n = itr->tmpvec[1]->n = itr->matches->n = itr->sub->n = 0;
	if (itr->start >= itr->len || itr->start < 0) return 0;
	while (itr->start < itr->len && itr->query[itr->start] > 3) ++itr->start; // skip ambiguous bases
	if (itr->start == itr->len) return 0;
	ori_start = itr->start;
	itr->start = bwt_smem1a(itr->bwt, itr->len, itr->query, ori_start, itr->min_intv, itr->max_intv, itr->matches, itr->tmpvec); // search for SMEM
	return itr->matches;
}

/***********************
 *** Extra functions ***
 ***********************/



