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

#include "bwa/bntseq.h"
#include "bwa/bwa.h"
#include "bwa/bwt.h"
#include "bwa/bwamem.h"
#include "bwa/kseq.h"
#include "bwa/ksort.h"
#include "bwa/kvec.h"
#include "bwa/utils.h"

#include "blaze/AccAgent.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.7.13-r1126-wrappered"
#endif
extern "C"{
KSEQ_DECLARE(gzFile)
}

extern blaze::AccAgent* agent;
const std::string acc_id = "SmithWaterman";

class ktp_aux_t{
public:
	kseq_t *ks, *ks2;
	mem_opt_t *opt;
	mem_pestat_t *pes0;
	int64_t n_processed;
	int copy_comment, actual_chunk_size;
	bwaidx_t *idx;
} ;


 class smem_aux_t {
 public:
    int id_read;
	bwtintv_v mem, mem1, *tmpv[2];
} ;

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

class mem_chain_v
{
public:
    size_t n, m; mem_chain_t *a;
 } ;

class MemChainVector
{
    public:
        int id_read;
        size_t n, m; mem_chain_t *a;
        MemChainVector()
        {
            id_read = 0;
            n = 0;
            m = 0;
            a = 0 ;
        }

};

struct mem_chainref_t {
  int64_t   rmax[2];
  uint8_t*  rseq;
  uint64_t* srt;
};

class MemAlignRegVector
{
    public:
        int id_read;
        size_t n, m;
        mem_alnreg_t *a;
};


 extern "C"{

int ksprintf(kstring_t *s, const char *fmt, ...);

void *kopen(const char *fn, int *_fd);

void mem_collect_intv(const mem_opt_t *opt, const bwt_t *bwt, int len, const uint8_t *seq, smem_aux_t *a);

mem_chain_v mem_chain(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, int len, const uint8_t *seq, void *buf);

int mem_chain_flt(const mem_opt_t *opt, int n_chn, mem_chain_t *a);

void mem_flt_chained_seeds(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, int n_chn, mem_chain_t *a);

void mem_chain2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, mem_alnreg_v *av);

int mem_sort_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, int n, mem_alnreg_t *a);

int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2], mem_alnreg_v a[2]);

int kclose(void *a);
}

int pre_process(int argc, char *argv[],ktp_aux_t *aux );

void reg_dump(mem_alnreg_v *alnreg,mem_alnreg_v *alnreg_hw,int batch_num);

void seq2intv(ktp_aux_t *aux,bseq1_t *seqs,smem_aux_t *SMEM);

mem_chain_v seq2chain(ktp_aux_t *aux, bseq1_t *seqs);

void chain2reg(ktp_aux_t *aux,bseq1_t *seqs,MemChainVector chn,mem_alnreg_v *alnreg);

void reg2sam(ktp_aux_t *aux,bseq1_t *seqs,int batch_num,int64_t n_processed,mem_alnreg_v *alnreg);

smem_aux_t *smem_aux_init();
void smem_aux_destroy(smem_aux_t *a);
int usage();

#endif // STRUCTURESNEW_H_INCLUDED
