#ifndef STRUCTURESNEW_H_INCLUDED
#define STRUCTURESNEW_H_INCLUDED
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
#include "bwt.h"
#include "bwamem.h"
#include "kvec.h"
#include "bntseq.h"
#include "kseq.h"
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.7.13-r1126-wrappered"
#endif
KSEQ_DECLARE(gzFile)





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
    int id_read;
    size_t n, m; mem_chain_t *a;
 } ;

int pre_process(int argc, char *argv[]);

int load_reads();

void seq2intv(bseq1_t *seqs,smem_aux_t *SMEM);

void intv2chain(smem_aux_t *SMEM,mem_chain_v *chain);

void chain2reg(mem_chain_v *chn,mem_alnreg_v *alnreg);

void reg2sam(mem_alnreg_v *alnreg);

void *kopen(const char *fn, int *_fd);
int kclose(void *a);
smem_aux_t *smem_aux_init();
void smem_aux_destroy(smem_aux_t *a);
int usage();
#endif // STRUCTURESNEW_H_INCLUDED
