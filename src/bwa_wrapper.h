#ifndef STRUCTURESNEW_H_INCLUDED
#define STRUCTURESNEW_H_INCLUDED

#include <zlib.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include <string>
#include <vector>

#include "bwa/bntseq.h"
#include "bwa/bwa.h"
#include "bwa/bwt.h"
#include "bwa/bwamem.h"
#include "bwa/kseq.h"
#include "bwa/ksort.h"
#include "bwa/ksw.h"
#include "bwa/kvec.h"
#include "bwa/utils.h"

#ifdef USE_HTSLIB
#include <htslib/sam.h>
#endif

extern "C"{
KSEQ_DECLARE(gzFile)
}

class ktp_aux_t{
public:
	kseq_t *ks, *ks2;
	mem_opt_t *opt;
	mem_pestat_t *pes0;
	int64_t n_processed;
	int copy_comment, actual_chunk_size;
	bwaidx_t *idx;
        bam_hdr_t *h;
#ifdef USE_HTSLIB
        samFile * out;
#endif
};

//class smem_aux_t {
// public:
//  int id_read;
//	bwtintv_v mem, mem1, *tmpv[2];
//};
typedef struct {
	bwtintv_v mem, mem1, *tmpv[2];
} smem_aux_t;

typedef struct {
  kstring_t name, comment, seq, qual;
  int last_char;
} kseq_new_t; 

/************
 * Chaining *
 ************/

typedef struct {
	int64_t rbeg;
	int32_t qbeg, len;
	int score;
} mem_seed_t; // unaligned memory

class mem_chain_t {
 public:
	int n, m, first, rid;
	uint32_t w:29, kept:2, is_alt:1;
	float frac_rep;
	int64_t pos;
	mem_seed_t *seeds;
} ;

class mem_chain_v {
 public:
  size_t n, m; mem_chain_t *a;
} ;

struct mem_chainref_t {
  int64_t   rmax[2];
  uint8_t*  rseq;
  uint64_t* srt;
};

extern "C"{

int ksprintf(kstring_t *s, const char *fmt, ...);

void *kopen(const char *fn, int *_fd);

void bwa_format_sam_hdr(const bntseq_t *bns, const char *rg_line, kstring_t *str);
bams_t *bams_init();
void bams_add(bams_t *bams, bam1_t *b);
void bams_destroy(bams_t *bams);

mem_chain_v mem_chain(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, int len, const uint8_t *seq, void *buf);

int mem_chain_flt(const mem_opt_t *opt, int n_chn, mem_chain_t *a);

void mem_flt_chained_seeds(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, int n_chn, mem_chain_t *a);

void mem_chain2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, mem_alnreg_v *av);

int mem_sort_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, int n, mem_alnreg_t *a);

int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2], mem_alnreg_v a[2] );

int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id);

void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m);

int kclose(void *a);
}

int pre_process(int argc, char *argv[], ktp_aux_t *aux, bool is_master);
int pack_bwa_mem_args(std::vector<const char*> & bwa_mem_args);

extern "C" {

void ks_introsort_mem_intv(size_t n, key_t array[]);
inline int get_rlen(int n_cigar, const uint32_t *cigar);
inline int mem_infer_dir(int64_t l_pac, int64_t b1, int64_t b2, int64_t *dist);
int mem_pair(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], bseq1_t s[2], mem_alnreg_v a[2], int id, int *sub, int *n_sub, int z[2], int n_pri[2]);
int mem_matesw(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], const mem_alnreg_t *a, int l_ms, const uint8_t *ms, mem_alnreg_v *ma);
int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a);
char **mem_gen_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_alnreg_v *a, int l_query, const char *query);

}

int kseq_read_new(kseq_new_t *seq_new, kseq_t *seq);

mem_chain_v mem_seq2chain_origin(ktp_aux_t *aux, bseq1_t *seqs);
mem_chain_v mem_seq2chain_wrapper(ktp_aux_t *aux, bseq1_t *seqs);

mem_chain_v mem_chain_new(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, int len, const uint8_t *seq, void *buf);

mem_chain_v mem_chain_postprocess(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, int len, const uint8_t *seq, void *buf);

void mem_collect_intv_new(const mem_opt_t *opt, const bwt_t *bwt, int len, const uint8_t *seq, smem_aux_t *a);

int bwt_smem1a_new(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_intv, uint64_t max_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2],int min_seed_len);

void mem_aln2sam_hts(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str, bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m_, bam_hdr_t *h);

void mem_reg2sam_hts(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m, bam_hdr_t *h);

int mem_sam_pe_hts(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2], mem_alnreg_v a[2], bam_hdr_t *header);

uint8_t *bns_fetch_seq_fpga(const bntseq_t *bns, const uint8_t *pac, int64_t *beg, int64_t mid, int64_t *end, int *rid);

void freeChains(mem_chain_v* chains, int batch_num);
void freeAligns(mem_alnreg_v* alnreg, int batch_num);
void freeSeqs(bseq1_t* seqs, int batch_num);

extern "C" {
smem_aux_t *smem_aux_init();
void smem_aux_destroy(smem_aux_t *a);
}
int usage();

#endif // STRUCTURESNEW_H_INCLUDED
