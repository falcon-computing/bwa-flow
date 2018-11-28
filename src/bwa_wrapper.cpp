#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif
extern "C" {
#include "bwa/kstring.h"
#include "bwa/bwamem.h"
#include "bwa/bntseq.h"
#include "bwa/ksw.h"
#include "bwa/kvec.h"
#include "bwa/ksort.h"
#include "bwa/kbtree.h"
#include "bwa/utils.h"
}
#include "bwa_wrapper.h"

#include "allocation_wrapper.h"


void bwa_format_sam_hdr(const bntseq_t *bns, const char *rg_line, kstring_t *str)
{
  int i;
  extern char *bwa_pg;
  str->l = 0; str->s = 0;
  ksprintf(str, "@HD\tVN:1.3\tSO:coordinate\n");
  for (i = 0; i < bns->n_seqs; ++i) 
    ksprintf(str, "@SQ\tSN:%s\tLN:%d\n", bns->anns[i].name, bns->anns[i].len);
  if (rg_line) ksprintf(str, "%s\n", rg_line);
  ksprintf(str, "%s\n", bwa_pg);
}

bams_t *bams_init() {
  return (bams_t *)calloc(1, sizeof(bams_t));
}

void bams_add(bams_t *bams, bam1_t *b) {
  if (bams->m == bams->l) {
    bams->m = bams->m ? bams->m << 1 : 4;
    bams->bams = (bam1_t **)realloc(bams->bams, sizeof(bam1_t*) * bams->m);
  }
  bams->bams[bams->l] = b;
  bams->l++;
}

void bams_destroy(bams_t *bams) {
  int i;
  for (i = 0; i < bams->l; i++) {
    bam_destroy1(bams->bams[i]);
  }
  free(bams->bams);
  free(bams);
}

extern "C" {
inline int ks_getc(kstream_t *ks);
inline int ks_getuntil(kstream_t *ks, int delimiter, kstring_t *str, int *dret);
int ks_getuntil_core(kstream_t *ks, int delimiter, kstring_t *str, int *dret, int append);
}
int kseq_read_new (kseq_new_t* seq_new, kseq_t *seq)
{
  int c;
  kstream_t *ks = seq->f;
  if (seq->last_char == 0) { /* then jump to the next header line */
    while ((c = ks_getc(ks)) != -1 && c != '>' && c != '@');
    if (c == -1) return -1; /* end of file */
    seq->last_char = c;
  } /* else: the first header char has been read in the previous call */
  seq_new->comment.l = seq_new->seq.l = seq_new->qual.l = 0; /* reset all members */
  if (ks_getuntil(ks, 0, &seq_new->name, &c) < 0) return -1; /* normal exit: EOF */
  if (c != '\n') ks_getuntil(ks, KS_SEP_LINE, &seq_new->comment, 0); /* read FASTA/Q comment */
  if (seq_new->seq.s == 0) { /* we can do this in the loop below, but that is slower */
    seq_new->seq.m = 256;
    seq_new->seq.s = (char*)malloc(seq_new->seq.m);
  }
  while ((c = ks_getc(ks)) != -1 && c != '>' && c != '+' && c != '@') {
    if (c == '\n') continue; /* skip empty lines */
    seq_new->seq.s[seq_new->seq.l++] = c; /* this is safe: we always have enough space for 1 char */
    ks_getuntil_core(ks, KS_SEP_LINE, &seq_new->seq, 0, 1); /* read the rest of the line */
  }
  if (c == '>' || c == '@') seq->last_char = c; /* the first header char has been read */
  if (seq_new->seq.l + 1 >= seq_new->seq.m) { /* seq_new->seq.s[seq_new->seq.l] below may be out of boundary */
    seq_new->seq.m = seq_new->seq.l + 2;
    kroundup32(seq_new->seq.m); /* rounded to the next closest 2^k */
    seq_new->seq.s = (char*)realloc(seq_new->seq.s, seq_new->seq.m);
  }
  seq_new->seq.s[seq_new->seq.l] = 0;     /* null terminated string */
  if (c != '+') return seq_new->seq.l; /* FASTA */
  if (seq_new->qual.m < seq_new->seq.m) { /* allocate memory for qual in case insufficient */
    seq_new->qual.m = seq_new->seq.m;
    seq_new->qual.s = (char*)realloc(seq_new->qual.s, seq_new->qual.m);
  }
  while ((c = ks_getc(ks)) != -1 && c != '\n'); /* skip the rest of '+' line */
  if (c == -1) return -2; /* error: no quality string */
  while (ks_getuntil_core(ks, KS_SEP_LINE, &seq_new->qual, 0, 1) >= 0 && seq_new->qual.l < seq_new->seq.l);
  seq->last_char = 0;     /* we have not come to the next header line */
  if (seq_new->seq.l != seq_new->qual.l) return -2; /* error: qual string is of a different length */
  return seq_new->seq.l;
}

mem_chain_v mem_seq2chain_origin(ktp_aux_t *aux, bseq1_t *seqs) {
  int i;
  mem_chain_v chn;
  for (i = 0; i < seqs->l_seq; ++i)
    seqs->seq[i] = seqs->seq[i] < 4? seqs->seq[i] : nst_nt4_table[(int)seqs->seq[i]];

  chn = mem_chain(aux->opt, aux->idx->bwt, aux->idx->bns, seqs->l_seq, (uint8_t*)seqs->seq, 0);
  chn.n = mem_chain_flt(aux->opt, chn.n, chn.a);
  mem_flt_chained_seeds(aux->opt, aux->idx->bns, aux->idx->pac, seqs->l_seq, (uint8_t*)seqs->seq, chn.n, chn.a);

  return chn;
}

mem_chain_v mem_seq2chain_wrapper(
    ktp_aux_t *aux,
    bseq1_t *seqs
) {
  int i;
  mem_chain_v chn;
  for (i = 0; i < seqs->l_seq; ++i) // convert to 2-bit encoding if we have not done so
    seqs->seq[i] = seqs->seq[i] < 4? seqs->seq[i] : nst_nt4_table[(int)seqs->seq[i]];

  chn = mem_chain_new(aux->opt, aux->idx->bwt, aux->idx->bns, seqs->l_seq, (uint8_t*)seqs->seq, 0);         // the 0 should be reconsidered
  chn.n = mem_chain_flt(aux->opt, chn.n, chn.a);
  mem_flt_chained_seeds(aux->opt, aux->idx->bns, aux->idx->pac, seqs->l_seq, (uint8_t*)seqs->seq, chn.n, chn.a);

  return chn;
}

typedef struct { kbnode_t *root; int off_key, off_ptr, ilen, elen; int n, t; int n_keys, n_nodes; } kbtree_chn_t;
extern "C" {
kbtree_chn_t *kb_init_chn(int size);
int test_and_merge(const mem_opt_t *opt, int64_t l_pac, mem_chain_t *c, const mem_seed_t *p, int seed_rid);
void kb_putp_chn(kbtree_chn_t *b, const mem_chain_t * __restrict k);
void kb_intervalp_chn(kbtree_chn_t *b, const mem_chain_t * __restrict k, mem_chain_t **lower, mem_chain_t **upper);
}
mem_chain_v mem_chain_new(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, int len, const uint8_t *seq, void *buf)
{
	int i, b, e, l_rep;
	int64_t l_pac = bns->l_pac;
	mem_chain_v chain;
	kbtree_t(chn) *tree;
	smem_aux_t *aux;

	kv_init(chain);
	if (len < opt->min_seed_len) return chain; // if the query is shorter than the seed length, no match
	tree = kb_init(chn, KB_DEFAULT_SIZE);

	aux = buf? (smem_aux_t*)buf : smem_aux_init();
	mem_collect_intv_new(opt, bwt, len, seq, aux);
	for (i = 0, b = e = l_rep = 0; i < aux->mem.n; ++i) { // compute frac_rep
		bwtintv_t *p = &aux->mem.a[i];
		int sb = (p->info>>32), se = (uint32_t)p->info;
		if (p->x[2] <= opt->max_occ) continue;
		if (sb > e) l_rep += e - b, b = sb, e = se;
		else e = e > se? e : se;
	}
	l_rep += e - b;
	for (i = 0; i < aux->mem.n; ++i) {
		bwtintv_t *p = &aux->mem.a[i];
		int step, count, slen = (uint32_t)p->info - (p->info>>32); // seed length
		int64_t k;
		// if (slen < opt->min_seed_len) continue; // ignore if too short or too repetitive
		step = p->x[2] > opt->max_occ? p->x[2] / opt->max_occ : 1;
		for (k = count = 0; k < p->x[2] && count < opt->max_occ; k += step, ++count) {
			mem_chain_t tmp, *lower, *upper;
			mem_seed_t s;
			int rid, to_add = 0;
			s.rbeg = tmp.pos = bwt_sa(bwt, p->x[0] + k); // this is the base coordinate in the forward-reverse reference
			s.qbeg = p->info>>32;
			s.score= s.len = slen;
			rid = bns_intv2rid(bns, s.rbeg, s.rbeg + s.len);
			if (rid < 0) continue; // bridging multiple reference sequences or the forward-reverse boundary; TODO: split the seed; don't discard it!!!
			if (kb_size(tree)) {
				kb_intervalp(chn, tree, &tmp, &lower, &upper); // find the closest chain
				if (!lower || !test_and_merge(opt, l_pac, lower, &s, rid)) to_add = 1;
			} else to_add = 1;
			if (to_add) { // add the seed as a new chain
				tmp.n = 1; tmp.m = 4;
				tmp.seeds = (mem_seed_t*)calloc(tmp.m, sizeof(mem_seed_t));
				tmp.seeds[0] = s;
				tmp.rid = rid;
				tmp.is_alt = !!bns->anns[rid].is_alt;
				kb_putp(chn, tree, &tmp);
			}
		}
	}
	if (buf == 0) smem_aux_destroy(aux);

	kv_resize(mem_chain_t, chain, kb_size(tree));

	#define traverse_func(p_) (chain.a[chain.n++] = *(p_))
	__kb_traverse(mem_chain_t, tree, traverse_func);
	#undef traverse_func

	for (i = 0; i < chain.n; ++i) chain.a[i].frac_rep = (float)l_rep / len;
	if (bwa_verbose >= 4) printf("* fraction of repetitive seeds: %.3f\n", (float)l_rep / len);

	kb_destroy(chn, tree);
	return chain;
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
  ks_introsort(mem_intv, a->mem.n, (int *)a->mem.a);
}

mem_chain_v mem_chain_postprocess(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, int len, const uint8_t *seq, void *buf)
{
#if 0
	mem_chain_v chain;
	kv_init(chain);
	if (len < opt->min_seed_len) return chain; // if the query is shorter than the seed length, no match

	smem_aux_t *aux;
	aux = buf? (smem_aux_t*)buf : smem_aux_init();
	//mem_collect_intv(opt, bwt, len, seq, aux);
	//use the new method
	mem_collect_intv_new(opt, bwt, len, seq, aux);
#endif
        smem_aux_t *a = (smem_aux_t *)buf;
        ks_introsort(mem_intv, a->mem.n, (int *)a->mem.a);
	int i, b, e, l_rep;
	int64_t l_pac = bns->l_pac;
	mem_chain_v chain;
	kbtree_t(chn) *tree;
	smem_aux_t *aux = (smem_aux_t*)buf;

	kv_init(chain);
	if (len < opt->min_seed_len) return chain; // if the query is shorter than the seed length, no match
	tree = kb_init(chn, KB_DEFAULT_SIZE);

	for (i = 0, b = e = l_rep = 0; i < aux->mem.n; ++i) { // compute frac_rep
		bwtintv_t *p = &aux->mem.a[i];
		int sb = (p->info>>32), se = (uint32_t)p->info;
		if (p->x[2] <= opt->max_occ) continue;
		if (sb > e) l_rep += e - b, b = sb, e = se;
		else e = e > se? e : se;
	}
	l_rep += e - b;
	for (i = 0; i < aux->mem.n; ++i) {
		bwtintv_t *p = &aux->mem.a[i];
		int step, count, slen = (uint32_t)p->info - (p->info>>32); // seed length
		int64_t k;
		// if (slen < opt->min_seed_len) continue; // ignore if too short or too repetitive
		step = p->x[2] > opt->max_occ? p->x[2] / opt->max_occ : 1;
		for (k = count = 0; k < p->x[2] && count < opt->max_occ; k += step, ++count) {
			mem_chain_t tmp, *lower, *upper;
			mem_seed_t s;
			int rid, to_add = 0;
			s.rbeg = tmp.pos = bwt_sa(bwt, p->x[0] + k); // this is the base coordinate in the forward-reverse reference
			s.qbeg = p->info>>32;
			s.score= s.len = slen;
			rid = bns_intv2rid(bns, s.rbeg, s.rbeg + s.len);
			if (rid < 0) continue; // bridging multiple reference sequences or the forward-reverse boundary; TODO: split the seed; don't discard it!!!
			if (kb_size(tree)) {
				kb_intervalp(chn, tree, &tmp, &lower, &upper); // find the closest chain
				if (!lower || !test_and_merge(opt, l_pac, lower, &s, rid)) to_add = 1;
			} else to_add = 1;
			if (to_add) { // add the seed as a new chain
				tmp.n = 1; tmp.m = 4;
				tmp.seeds = (mem_seed_t*)calloc(tmp.m, sizeof(mem_seed_t));
				tmp.seeds[0] = s;
				tmp.rid = rid;
				tmp.is_alt = !!bns->anns[rid].is_alt;
				kb_putp(chn, tree, &tmp);
			}
		}
	}
	if (buf == 0) smem_aux_destroy(aux);

	kv_resize(mem_chain_t, chain, kb_size(tree));

	#define traverse_func(p_) (chain.a[chain.n++] = *(p_))
	__kb_traverse(mem_chain_t, tree, traverse_func);
	#undef traverse_func

	for (i = 0; i < chain.n; ++i) chain.a[i].frac_rep = (float)l_rep / len;
	if (bwa_verbose >= 4) printf("* fraction of repetitive seeds: %.3f\n", (float)l_rep / len);

	kb_destroy(chn, tree);
	return chain;
}

int bwt_smem1a_new(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_intv, uint64_t max_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2],int min_seed_len)
{
  int i, j, c, ret;
  bwtintv_t ik, ok[4];
  bwtintv_v a[2], *curr;
  int start = x;
  int stop = x;
  bwtintv_v *back_intv; 
  bwtintv_t temp;
  int k = 0;
  int m = 0;
  int isbreak = 0;
  // forward search
  if (q[x] > 3) return x + 1;
  if (min_intv < 1) min_intv = 1; // the interval size should be at least 1
  kv_init(a[0]); kv_init(a[1]);
  back_intv = &a[0];
  curr = &a[1];
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
            isbreak = 1;
          }
        }
        if(isbreak == 1){
          isbreak = 0;
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
    if (i < curr->n) {
      int max_len =(temp.info>>32) + (int)curr->a[i].info;
      while (max_len < min_seed_len && i < curr->n){
        i = i + 1;
        if (i < curr->n)
          stop = (int)curr->a[i].info;
        max_len = (temp.info>>32) + stop;
      }
    }
    if(i >= curr->n && (int)temp.info - (temp.info>>32) >= min_seed_len )
      kv_push(bwtintv_t,*mem,temp);
  }
  free (back_intv->a);
  free (a[1].a);
  //if (tmpvec == 0 || tmpvec[0] == 0) free(a[0].a);
  //if (tmpvec == 0 || tmpvec[1] == 0) free(a[1].a);
  return ret;
}

void mem_aln2sam_hts(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str, bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m_, bam_hdr_t *h)
{
	int i, l_name;
	mem_aln_t ptmp = list[which], *p = &ptmp, mtmp, *m = 0; // make a copy of the alignment to convert

	if (m_) mtmp = *m_, m = &mtmp;
	// set flag
	p->flag |= m? 0x1 : 0; // is paired in sequencing
	p->flag |= p->rid < 0? 0x4 : 0; // is mapped
	p->flag |= m && m->rid < 0? 0x8 : 0; // is mate mapped
	if (p->rid < 0 && m && m->rid >= 0) // copy mate to alignment
		p->rid = m->rid, p->pos = m->pos, p->is_rev = m->is_rev, p->n_cigar = 0;
	if (m && m->rid < 0 && p->rid >= 0) // copy alignment to mate
		m->rid = p->rid, m->pos = p->pos, m->is_rev = p->is_rev, m->n_cigar = 0;
	p->flag |= p->is_rev? 0x10 : 0; // is on the reverse strand
	p->flag |= m && m->is_rev? 0x20 : 0; // is mate on the reverse strand

	// print up to CIGAR
	l_name = strlen(s->name);
	ks_resize(str, str->l + s->l_seq + l_name + (s->qual? s->l_seq : 0) + 20);
	kputsn(s->name, l_name, str); kputc('\t', str); // QNAME
	kputw((p->flag&0xffff) | (p->flag&0x10000? 0x100 : 0), str); kputc('\t', str); // FLAG
	if (p->rid >= 0) { // with coordinate
		kputs(bns->anns[p->rid].name, str); kputc('\t', str); // RNAME
		kputl(p->pos + 1, str); kputc('\t', str); // POS
		kputw(p->mapq, str); kputc('\t', str); // MAPQ
		if (p->n_cigar) { // aligned
			for (i = 0; i < p->n_cigar; ++i) {
				int c = p->cigar[i]&0xf;
				if (!(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt && (c == 3 || c == 4))
					c = which? 4 : 3; // use hard clipping for supplementary alignments
				kputw(p->cigar[i]>>4, str); kputc("MIDSH"[c], str);
			}
		} else kputc('*', str); // having a coordinate but unaligned (e.g. when copy_mate is true)
	} else kputsn("*\t0\t0\t*", 7, str); // without coordinte
	kputc('\t', str);

	// print the mate position if applicable
	if (m && m->rid >= 0) {
		if (p->rid == m->rid) kputc('=', str);
		else kputs(bns->anns[m->rid].name, str);
		kputc('\t', str);
		kputl(m->pos + 1, str); kputc('\t', str);
		if (p->rid == m->rid) {
			int64_t p0 = p->pos + (p->is_rev? get_rlen(p->n_cigar, p->cigar) - 1 : 0);
			int64_t p1 = m->pos + (m->is_rev? get_rlen(m->n_cigar, m->cigar) - 1 : 0);
			if (m->n_cigar == 0 || p->n_cigar == 0) kputc('0', str);
			else kputl(-(p0 - p1 + (p0 > p1? 1 : p0 < p1? -1 : 0)), str);
		} else kputc('0', str);
	} else kputsn("*\t0\t0", 5, str);
	kputc('\t', str);

	// print SEQ and QUAL
	if (p->flag & 0x100) { // for secondary alignments, don't write SEQ and QUAL
		kputsn("*\t*", 3, str);
	} else if (!p->is_rev) { // the forward strand
		int i, qb = 0, qe = s->l_seq;
		if (p->n_cigar && which && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) { // have cigar && not the primary alignment && not softclip all
			if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qb += p->cigar[0]>>4;
			if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qe -= p->cigar[p->n_cigar-1]>>4;
		}
		ks_resize(str, str->l + (qe - qb) + 1);
		for (i = qb; i < qe; ++i) str->s[str->l++] = "ACGTN"[(int)s->seq[i]];
		kputc('\t', str);
		if (s->qual) { // printf qual
			ks_resize(str, str->l + (qe - qb) + 1);
			for (i = qb; i < qe; ++i) str->s[str->l++] = s->qual[i];
			str->s[str->l] = 0;
		} else kputc('*', str);
	} else { // the reverse strand
		int i, qb = 0, qe = s->l_seq;
		if (p->n_cigar && which && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) {
			if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qe -= p->cigar[0]>>4;
			if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qb += p->cigar[p->n_cigar-1]>>4;
		}
		ks_resize(str, str->l + (qe - qb) + 1);
		for (i = qe-1; i >= qb; --i) str->s[str->l++] = "TGCAN"[(int)s->seq[i]];
		kputc('\t', str);
		if (s->qual) { // printf qual
			ks_resize(str, str->l + (qe - qb) + 1);
			for (i = qe-1; i >= qb; --i) str->s[str->l++] = s->qual[i];
			str->s[str->l] = 0;
		} else kputc('*', str);
	}

	// print optional tags
	if (p->n_cigar) {
		kputsn("\tNM:i:", 6, str); kputw(p->NM, str);
		kputsn("\tMD:Z:", 6, str); kputs((char*)(p->cigar + p->n_cigar), str);
	}
	if (p->score >= 0) { kputsn("\tAS:i:", 6, str); kputw(p->score, str); }
	if (p->sub >= 0) { kputsn("\tXS:i:", 6, str); kputw(p->sub, str); }
	if (bwa_rg_id[0]) { kputsn("\tRG:Z:", 6, str); kputs(bwa_rg_id, str); }
	if (!(p->flag & 0x100)) { // not multi-hit
		for (i = 0; i < n; ++i)
			if (i != which && !(list[i].flag&0x100)) break;
		if (i < n) { // there are other primary hits; output them
			kputsn("\tSA:Z:", 6, str);
			for (i = 0; i < n; ++i) {
				const mem_aln_t *r = &list[i];
				int k;
				if (i == which || (r->flag&0x100)) continue; // proceed if: 1) different from the current; 2) not shadowed multi hit
				kputs(bns->anns[r->rid].name, str); kputc(',', str);
				kputl(r->pos+1, str); kputc(',', str);
				kputc("+-"[r->is_rev], str); kputc(',', str);
				for (k = 0; k < r->n_cigar; ++k) {
					kputw(r->cigar[k]>>4, str); kputc("MIDSH"[r->cigar[k]&0xf], str);
				}
				kputc(',', str); kputw(r->mapq, str);
				kputc(',', str); kputw(r->NM, str);
				kputc(';', str);
			}
		}
		if (p->alt_sc > 0)
			ksprintf(str, "\tpa:f:%.3f", (double)p->score / p->alt_sc);
	}
	if (p->XA) { kputsn("\tXA:Z:", 6, str); kputs(p->XA, str); }
	if (s->comment) { kputc('\t', str); kputs(s->comment, str); }
	if ((opt->flag&MEM_F_REF_HDR) && p->rid >= 0 && bns->anns[p->rid].anno != 0 && bns->anns[p->rid].anno[0] != 0) {
		int tmp;
		kputsn("\tXR:Z:", 6, str);
		tmp = str->l;
		kputs(bns->anns[p->rid].anno, str);
		for (i = tmp; i < str->l; ++i) // replace TAB in the comment to SPACE
			if (str->s[i] == '\t') str->s[i] = ' ';
	}
  
#ifdef USE_HTSLIB
  bam1_t *b = bam_init1();
  if (sam_parse1(str, h, b) < 0) {
    fprintf(stderr, "sam_parse1() error!\n");
    // TODO(mhhuang): needs error handling here and exit!
  }
  bams_add(s->bams, b);
  str->l = 0; 
  //str->s = 0;
#else
	kputc('\n', str);
#endif
}

void mem_reg2sam_hts(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m, bam_hdr_t *h)
{
	kstring_t str;
	kvec_t(mem_aln_t) aa;
	int k, l;
	char **XA = 0;

	if (!(opt->flag & MEM_F_ALL))
		XA = mem_gen_alt(opt, bns, pac, a, s->l_seq, s->seq);
	kv_init(aa);
	str.l = str.m = 0; str.s = 0;
	for (k = l = 0; k < a->n; ++k) {
		mem_alnreg_t *p = &a->a[k];
		mem_aln_t *q;
		if (p->score < opt->T) continue;
		if (p->secondary >= 0 && (p->is_alt || !(opt->flag&MEM_F_ALL))) continue;
		if (p->secondary >= 0 && p->secondary < INT_MAX && p->score < a->a[p->secondary].score * opt->drop_ratio) continue;
		q = kv_pushp(mem_aln_t, aa);
		*q = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, p);
		assert(q->rid >= 0); // this should not happen with the new code
		q->XA = XA? XA[k] : 0;
		q->flag |= extra_flag; // flag secondary
		if (p->secondary >= 0) q->sub = -1; // don't output sub-optimal score
		if (l && p->secondary < 0) // if supplementary
			q->flag |= (opt->flag&MEM_F_NO_MULTI)? 0x10000 : 0x800;
		if (l && !p->is_alt && q->mapq > aa.a[0].mapq) q->mapq = aa.a[0].mapq;
		++l;
	}
	if (aa.n == 0) { // no alignments good enough; then write an unaligned record
		mem_aln_t t;
		t = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, 0);
		t.flag |= extra_flag;
		mem_aln2sam_hts(opt, bns, &str, s, 1, &t, 0, m, h);
	} else {
		for (k = 0; k < aa.n; ++k)
			mem_aln2sam_hts(opt, bns, &str, s, aa.n, aa.a, k, m, h);
		for (k = 0; k < aa.n; ++k) free(aa.a[k].cigar);
		free(aa.a);
	}
#ifdef USE_HTSLIB
  free(str.s);
#else
	s->sam = str.s;
#endif
	if (XA) {
		for (k = 0; k < a->n; ++k) free(XA[k]);
		free(XA);
	}
}

#define raw_mapq(diff, a) ((int)(6.02 * (diff) / (a) + .499))
int mem_sam_pe_hts(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2], mem_alnreg_v a[2], bam_hdr_t *header)
{
	//extern int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id); // declared in "bwa_wrapper.h"
	//extern void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m, bam_hdr_t *h); // not called here

	int n = 0, i, j, z[2], o, subo, n_sub, extra_flag = 1, n_pri[2], n_aa[2];
	kstring_t str;
	mem_aln_t h[2], g[2], aa[2][2];

	str.l = str.m = 0; str.s = 0;
	memset(h, 0, sizeof(mem_aln_t) * 2);
	memset(g, 0, sizeof(mem_aln_t) * 2);
	n_aa[0] = n_aa[1] = 0;
	if (!(opt->flag & MEM_F_NO_RESCUE)) { // then perform SW for the best alignment
		mem_alnreg_v b[2];
		kv_init(b[0]); kv_init(b[1]);
		for (i = 0; i < 2; ++i)
			for (j = 0; j < a[i].n; ++j)
				if (a[i].a[j].score >= a[i].a[0].score  - opt->pen_unpaired)
					kv_push(mem_alnreg_t, b[i], a[i].a[j]);
		for (i = 0; i < 2; ++i)
			for (j = 0; j < b[i].n && j < opt->max_matesw; ++j)
				n += mem_matesw(opt, bns, pac, pes, &b[i].a[j], s[!i].l_seq, (uint8_t*)s[!i].seq, &a[!i]);
		free(b[0].a); free(b[1].a);
	}
	n_pri[0] = mem_mark_primary_se(opt, a[0].n, a[0].a, id<<1|0);
	n_pri[1] = mem_mark_primary_se(opt, a[1].n, a[1].a, id<<1|1);
	if (opt->flag&MEM_F_NOPAIRING) goto no_pairing;
	// pairing single-end hits
	if (n_pri[0] && n_pri[1] && (o = mem_pair(opt, bns, pac, pes, s, a, id, &subo, &n_sub, z, n_pri)) > 0) {
		int is_multi[2], q_pe, score_un, q_se[2];
		char **XA[2];
		// check if an end has multiple hits even after mate-SW
		for (i = 0; i < 2; ++i) {
			for (j = 1; j < n_pri[i]; ++j)
				if (a[i].a[j].secondary < 0 && a[i].a[j].score >= opt->T) break;
			is_multi[i] = j < n_pri[i]? 1 : 0;
		}
		if (is_multi[0] || is_multi[1]) goto no_pairing; // TODO: in rare cases, the true hit may be long but with low score
		// compute mapQ for the best SE hit
		score_un = a[0].a[0].score + a[1].a[0].score - opt->pen_unpaired;
		//q_pe = o && subo < o? (int)(MEM_MAPQ_COEF * (1. - (double)subo / o) * log(a[0].a[z[0]].seedcov + a[1].a[z[1]].seedcov) + .499) : 0;
		subo = subo > score_un? subo : score_un;
		q_pe = raw_mapq(o - subo, opt->a);
		if (n_sub > 0) q_pe -= (int)(4.343 * log(n_sub+1) + .499);
		if (q_pe < 0) q_pe = 0;
		if (q_pe > 60) q_pe = 60;
		q_pe = (int)(q_pe * (1. - .5 * (a[0].a[0].frac_rep + a[1].a[0].frac_rep)) + .499);
		// the following assumes no split hits
		if (o > score_un) { // paired alignment is preferred
			mem_alnreg_t *c[2];
			c[0] = &a[0].a[z[0]]; c[1] = &a[1].a[z[1]];
			for (i = 0; i < 2; ++i) {
				if (c[i]->secondary >= 0)
					c[i]->sub = a[i].a[c[i]->secondary].score, c[i]->secondary = -2;
				q_se[i] = mem_approx_mapq_se(opt, c[i]);
			}
			q_se[0] = q_se[0] > q_pe? q_se[0] : q_pe < q_se[0] + 40? q_pe : q_se[0] + 40;
			q_se[1] = q_se[1] > q_pe? q_se[1] : q_pe < q_se[1] + 40? q_pe : q_se[1] + 40;
			extra_flag |= 2;
			// cap at the tandem repeat score
			q_se[0] = q_se[0] < raw_mapq(c[0]->score - c[0]->csub, opt->a)? q_se[0] : raw_mapq(c[0]->score - c[0]->csub, opt->a);
			q_se[1] = q_se[1] < raw_mapq(c[1]->score - c[1]->csub, opt->a)? q_se[1] : raw_mapq(c[1]->score - c[1]->csub, opt->a);
		} else { // the unpaired alignment is preferred
			z[0] = z[1] = 0;
			q_se[0] = mem_approx_mapq_se(opt, &a[0].a[0]);
			q_se[1] = mem_approx_mapq_se(opt, &a[1].a[0]);
		}
		for (i = 0; i < 2; ++i) {
			int k = a[i].a[z[i]].secondary_all;
			if (k >= 0 && k < n_pri[i]) { // switch secondary and primary if both of them are non-ALT
				assert(a[i].a[k].secondary_all < 0);
				for (j = 0; j < a[i].n; ++j)
					if (a[i].a[j].secondary_all == k || j == k)
						a[i].a[j].secondary_all = z[i];
				a[i].a[z[i]].secondary_all = -1;
			}
		}
		if (!(opt->flag & MEM_F_ALL)) {
			for (i = 0; i < 2; ++i)
				XA[i] = mem_gen_alt(opt, bns, pac, &a[i], s[i].l_seq, s[i].seq);
		} else XA[0] = XA[1] = 0;
		// write SAM
		for (i = 0; i < 2; ++i) {
			h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, &a[i].a[z[i]]);
			h[i].mapq = q_se[i];
			h[i].flag |= 0x40<<i | extra_flag;
			h[i].XA = XA[i]? XA[i][z[i]] : 0;
			aa[i][n_aa[i]++] = h[i];
			if (n_pri[i] < a[i].n) { // the read has ALT hits
				mem_alnreg_t *p = &a[i].a[n_pri[i]];
				if (p->score < opt->T || p->secondary >= 0 || !p->is_alt) continue;
				g[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, p);
				g[i].flag |= 0x800 | 0x40<<i | extra_flag;
				g[i].XA = XA[i]? XA[i][n_pri[i]] : 0;
				aa[i][n_aa[i]++] = g[i];
			}
		}
		for (i = 0; i < n_aa[0]; ++i)
			mem_aln2sam_hts(opt, bns, &str, &s[0], n_aa[0], aa[0], i, &h[1], header); // write read1 hits
#ifndef USE_HTSLIB
		s[0].sam = strdup(str.s);
#endif
    str.l = 0;
		for (i = 0; i < n_aa[1]; ++i)
			mem_aln2sam_hts(opt, bns, &str, &s[1], n_aa[1], aa[1], i, &h[0], header); // write read2 hits
#ifndef USE_HTSLIB
		s[1].sam = str.s;
#else
    free(str.s);
#endif
		if (strcmp(s[0].name, s[1].name) != 0) err_fatal(__func__, "paired reads have different names: \"%s\", \"%s\"\n", s[0].name, s[1].name);
		// free
		for (i = 0; i < 2; ++i) {
			free(h[i].cigar); free(g[i].cigar);
			if (XA[i] == 0) continue;
			for (j = 0; j < a[i].n; ++j) free(XA[i][j]);
			free(XA[i]);
		}
	} else goto no_pairing;
	return n;

no_pairing:
	for (i = 0; i < 2; ++i) {
		int which = -1;
		if (a[i].n) {
			if (a[i].a[0].score >= opt->T) which = 0;
			else if (n_pri[i] < a[i].n && a[i].a[n_pri[i]].score >= opt->T)
				which = n_pri[i];
		}
		if (which >= 0) h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, &a[i].a[which]);
		else h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, 0);
	}
	if (!(opt->flag & MEM_F_NOPAIRING) && h[0].rid == h[1].rid && h[0].rid >= 0) { // if the top hits from the two ends constitute a proper pair, flag it.
		int64_t dist;
		int d;
		d = mem_infer_dir(bns->l_pac, a[0].a[0].rb, a[1].a[0].rb, &dist);
		if (!pes[d].failed && dist >= pes[d].low && dist <= pes[d].high) extra_flag |= 2;
	}
	mem_reg2sam_hts(opt, bns, pac, &s[0], &a[0], 0x41|extra_flag, &h[1], header);
	mem_reg2sam_hts(opt, bns, pac, &s[1], &a[1], 0x81|extra_flag, &h[0], header);
	if (strcmp(s[0].name, s[1].name) != 0) err_fatal(__func__, "paired reads have different names: \"%s\", \"%s\"\n", s[0].name, s[1].name);
	free(h[0].cigar); free(h[1].cigar);
	return n;
}

uint8_t *bns_fetch_seq_fpga(const bntseq_t *bns, const uint8_t *pac, int64_t *beg, int64_t mid, int64_t *end, int *rid)
{
	int64_t far_beg, far_end, len;
	int is_rev;
	uint8_t *seq;

	if (*end < *beg) *end ^= *beg, *beg ^= *end, *end ^= *beg; // if end is smaller, swap
	assert(*beg <= mid && mid < *end);
	*rid = bns_pos2rid(bns, bns_depos(bns, mid, &is_rev));
	far_beg = bns->anns[*rid].offset;
	far_end = far_beg + bns->anns[*rid].len;
	if (is_rev) { // flip to the reverse strand
		int64_t tmp = far_beg;
		far_beg = (bns->l_pac<<1) - far_end;
		far_end = (bns->l_pac<<1) - tmp;
	}
	*beg = *beg > far_beg? *beg : far_beg;
	*end = *end < far_end? *end : far_end;
	return seq;
}

void freeChains(mem_chain_v* chains, int batch_num) {
  // Free the chains
  for (int seq_id = 0; seq_id < batch_num; seq_id++) {
    for (int chain_id = 0; chain_id < chains[seq_id].n; chain_id++) {
      if (chains[seq_id].a[chain_id].seeds && chains[seq_id].a[chain_id].n > 0)
        free(chains[seq_id].a[chain_id].seeds);
    }
    if (chains[seq_id].a && chains[seq_id].n > 0)
      free(chains[seq_id].a);
  }
  free(chains);
}

void freeAligns(mem_alnreg_v* alnreg, int batch_num) {
  for (int seq_id = 0; seq_id < batch_num; seq_id++) {
    if (alnreg[seq_id].a && alnreg[seq_id].n > 0)
      free(alnreg[seq_id].a); 
  }
  free(alnreg);
}

void freeSeqs(bseq1_t* seqs, int batch_num) {
  for (int seq_id = 0; seq_id < batch_num; seq_id++) {
    free(seqs[seq_id].name); 
    if (seqs[seq_id].comment) free(seqs[seq_id].comment);
    free(seqs[seq_id].seq); 
    free(seqs[seq_id].qual); 
#ifdef USE_HTSLIB
    if (seqs[seq_id].bams) free(seqs[seq_id].bams);
#else
    if (seqs[seq_id].sam) free(seqs[seq_id].sam);
#endif
  }
  free(seqs);
}
