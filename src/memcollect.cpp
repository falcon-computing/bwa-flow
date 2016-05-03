#include <stdlib.h>
#include "bwa_wrapper.h"
#include "bwa/bwamem.h"
#include "bwa/bwt.h"
#include "bwa/ksort.h"
#include "kvec.h"
smem_aux_t *smem_aux_init()
{
	smem_aux_t *a;
	a = (smem_aux_t *)calloc(1, sizeof(smem_aux_t));
	a->tmpv[0] = (bwtintv_v*)calloc(1, sizeof(bwtintv_v));
	a->tmpv[1] = (bwtintv_v*)calloc(1, sizeof(bwtintv_v));
	return a;
}

#define intv_lt(a, b) ((a).info < (b).info)
KSORT_INIT(mem_intv, bwtintv_t, intv_lt)

// NOTE: $max_intv is not currently used in BWA-MEM
int bwt_smem1a_new(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_intv, uint64_t max_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2],int min_seed_len)
{
	int i, j, c, ret;
	bwtintv_t ik, ok[4];
	bwtintv_v a[2], *curr;
        int start = x;
        int stop = x;
        bwtintv_v *back_intv = new bwtintv_v;
        kv_init(*back_intv);
        bwtintv_t temp;
        int k = 0;
        int m = 0;
        bool isbreak = false;
        // forward search
	if (q[x] > 3) return x + 1;
	if (min_intv < 1) min_intv = 1; // the interval size should be at least 1
	kv_init(a[0]); kv_init(a[1]);
	curr = tmpvec && tmpvec[1]? tmpvec[1] : &a[1];
	bwt_set_intv(bwt, q[x], ik); // the initial interval of a single base
	ik.info = x + 1;

	for (i = x + 1, curr->n = 0; i < len; ++i) { // forward search
		if (ik.x[2] < max_intv) { // an interval small enough
			kv_push(bwtintv_t, *curr, ik);
			break;
		} else if (q[i] < 4) { // an A/C/G/T base
			c = 3 - q[i]; // complement of q[i]
			bwt_extend(bwt, &ik, ok, 0);
			if (ok[c].x[2] != ik.x[2]) { // change of the interval size
				kv_push(bwtintv_t, *curr, ik);
				if (ok[c].x[2] < min_intv) break; // the interval size is too small to be extended further
			}
			ik = ok[c]; ik.info = i + 1;
		} else { // an ambiguous base
			kv_push(bwtintv_t, *curr, ik);
			break; // always terminate extension at an ambiguous base; in this case, i<len always stands
		}
	}
	if (i == len) kv_push(bwtintv_t, *curr, ik); // push the last interval if we reach the end
	ret = curr->a[curr->n-1].info; // this will be the returned value
        // the new backward search 
        i = 0;
        while (i < curr->n){
          ik = curr->a[i];
          ik.info |=(uint64_t)(x)<<32;
          // backenlarge
          if(back_intv->n ==0 || stop - start >= 3){
             back_intv->n =0;
             kv_push(bwtintv_t, *back_intv,ik);
             for (k = x-1;k >= -1;--k){
               c = k<0?-1: q[k] < 4? q[k] : -1;
               if(c >=0 && ik.x[2] >= max_intv){
                 bwt_extend(bwt,&ik,ok,1);
                 if(ok[c].x[2] < min_intv) break;
                 ik = ok[c];
                 ik.info =curr->a[i].info|(uint64_t)(k)<<32;
                 kv_push(bwtintv_t, *back_intv,ik);
               }
               else {
                 break;
               }
             }
             start = (int)curr->a[i].info;
             if(i == curr->n - 1) stop = len;
             else{
               stop = (int)curr->a[i+1].info;
             }
             if(i==0)
                temp = ik;
             else if(ik.info>>32 > temp.info>>32 &&                             
                     (int)temp.info - (temp.info>>32) >= min_seed_len ){
                kv_push(bwtintv_t,*mem,temp);
                temp = ik;
             }
             else 
                temp = ik;
          }
          //forwardenlarge 
          else{
            stop = (int)curr->a[i].info;
            for (k = back_intv->n - 1; k >=0; k--){
              ik = back_intv->a[k];
              for (m = start + 1; m <= stop; m++){
                c =3 - q[m-1];
                bwt_extend(bwt,&ik,ok,0);
                if(ok[c].x[2] < min_intv)
                  break;
                ik = ok[c];
                if(m == stop){
                  ik.info = curr->a[i].info | (uint64_t)(x-k)<<32;
                  isbreak = true;
                }
              }
              if(isbreak == true){
                isbreak = false;
                if((x-k) > temp.info>>32 &&
                   (int)temp.info - (temp.info>>32) >= min_seed_len ){
                  kv_push(bwtintv_t,*mem,temp);
                  temp = ik;
                } 
                else
                  temp = ik;
                break;
              }
            }
          } 
          i = i+1; 
          // compute the max len of next iteration
          int max_len =(temp.info>>32) + (int)curr->a[i].info;
          while (max_len < min_seed_len && i < curr->n){
            i = i + 1;
            stop = (int)curr->a[i].info;
            max_len = (temp.info>>32) + stop;
          }
          if(i >= curr->n && (int)temp.info - (temp.info>>32) >= min_seed_len )
            kv_push(bwtintv_t,*mem,temp);
        }
        free (back_intv->a);
	if (tmpvec == 0 || tmpvec[0] == 0) free(a[0].a);
	if (tmpvec == 0 || tmpvec[1] == 0) free(a[1].a);
	return ret;
}



 void mem_collect_intv_new(const mem_opt_t *opt, const bwt_t *bwt, int len, const uint8_t *seq, smem_aux_t *a)
{
	int i, k, x = 0, old_n;
	int start_width = 1;
	int split_len = (int)(opt->min_seed_len * opt->split_factor + .499);
	a->mem.n = 0;
	// first pass: find all SMEMs
	while (x < len) {
		if (seq[x] < 4) {
			x = bwt_smem1a_new(bwt, len, seq, x, start_width,0, &a->mem, a->tmpv, opt->min_seed_len);
		} else ++x;
	}
	// second pass: find MEMs inside a long SMEM
	old_n = a->mem.n;
	for (k = 0; k < old_n; ++k) {
		bwtintv_t *p = &a->mem.a[k];
		int start = p->info>>32, end = (int32_t)p->info;
		if (end - start < split_len || p->x[2] > opt->split_width) continue;
		bwt_smem1a_new(bwt, len, seq, (start + end)>>1, p->x[2]+1,0, &a->mem, a->tmpv, opt->min_seed_len);
	}
	// third pass: LAST-like
	if (opt->max_mem_intv > 0) {
		x = 0;
		while (x < len) {
			if (seq[x] < 4) {
				if (1) {
					bwtintv_t m;
					x = bwt_seed_strategy1(bwt, len, seq, x, opt->min_seed_len, opt->max_mem_intv, &m);
					if (m.x[2] > 0) kv_push(bwtintv_t, a->mem, m);
				} else { // for now, we never come to this block which is slower
					x = bwt_smem1a(bwt, len, seq, x, start_width, opt->max_mem_intv, &a->mem1, a->tmpv);
					for (i = 0; i < a->mem1.n; ++i)
						kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
				}
			} else ++x;
		}
	}
	// sort
	ks_introsort(mem_intv, a->mem.n, a->mem.a);


}



