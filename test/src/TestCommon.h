#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include <unistd.h>
#include <dlfcn.h>

#include <boost/filesystem.hpp>
#include <cstdint>
#include <string>
#include <stdexcept>

#include <glog/logging.h>
#include <gtest/gtest.h>

#include "bwa/bntseq.h"
#include "bwa/bwa.h"
#include "bwa/bwamem.h"
#include "bwa/kseq.h"

#include "bwa_wrapper.h"

#ifdef USE_MPI
#include "mpi.h"
#include "MPIChannel.h"
#endif

// Global variables
extern char *bwa_pg;
extern gzFile fp_idx, fp2_read2;
extern void *ko_read1, *ko_read2;
extern ktp_aux_t *aux;

extern bseq1_t* g_seqs;
extern int g_batch_num;

static void dup_bseq1_t(bseq1_t* out, bseq1_t* in, int size) {
  for (int i = 0; i < size; i++) {
    out[i].name    = strdup(in[i].name);
    if (in[i].comment) {
      out[i].comment = strdup(in[i].comment);
    }
    else {
      out[i].comment = NULL;
    }
    out[i].seq     = strdup(in[i].seq);
    out[i].qual    = strdup(in[i].qual);
    out[i].l_seq   = strlen(in[i].seq);
  }
}

static bseq1_t* bwa_mem(bseq1_t* input, int batch_num) {

  // output bseq containing bams
  bseq1_t* seqs = (bseq1_t*)malloc(batch_num*sizeof(bseq1_t));
  dup_bseq1_t(seqs, input, batch_num);

  // temporary regions
  mem_alnreg_v* alnreg = new mem_alnreg_v[batch_num];

  for (int i = 0; i < batch_num; i++) {
    mem_chain_v chains = mem_seq2chain_wrapper(aux, &seqs[i]);
    kv_init(alnreg[i]);
    for (int j = 0; j < chains.n; j++) {
      mem_chain2aln(
          aux->opt, 
          aux->idx->bns, 
          aux->idx->pac,
          seqs[i].l_seq,
          (uint8_t*)seqs[i].seq,
          &chains.a[j],
          &alnreg[i]);

      free(chains.a[j].seeds);
    }
    free(chains.a);

    // Post-process each chain before output
    alnreg[i].n = mem_sort_dedup_patch(
        aux->opt, 
        aux->idx->bns, 
        aux->idx->pac, 
        (uint8_t*)seqs[i].seq, 
        alnreg[i].n, 
        alnreg[i].a);

    for (int j = 0; j < alnreg[i].n; j++) {
      mem_alnreg_t *p = &alnreg[i].a[j];
      if (p->rid >= 0 && aux->idx->bns->anns[p->rid].is_alt)
        p->is_alt = 1;
    }
  }

  // pair-end data
  if(aux->opt->flag&MEM_F_PE) {
    mem_pestat_t pes[4];
    mem_pestat(aux->opt, aux->idx->bns->l_pac, batch_num, alnreg, pes);

    for (int i = 0; i < batch_num/2; i++) {
      seqs[i<<1].bams = bams_init();
      seqs[1+(i<<1)].bams = bams_init();
      mem_sam_pe_hts(
          aux->opt,
          aux->idx->bns,
          aux->idx->pac,
          pes, i,
          &seqs[i<<1],
          &alnreg[i<<1],
          aux->h);
    }
  }
  else { // single-end data
    for (int i = 0; i < batch_num; i++) {
      seqs[i].bams = bams_init();
      mem_mark_primary_se(
          aux->opt,
          alnreg[i].n,
          alnreg[i].a,
          i);
      mem_reg2sam_hts(
          aux->opt,
          aux->idx->bns,
          aux->idx->pac,
          &seqs[i],
          &alnreg[i],
          0, 0, aux->h);
    }
  }
  freeAligns(alnreg, batch_num);
  for (int i = 0; i < batch_num; i++) {
    free(seqs[i].name);
    free(seqs[i].comment);
    free(seqs[i].seq);
    free(seqs[i].qual);
  }
  return seqs;
}



class BaseTests : public ::testing::Test {
  protected:
    virtual void SetUp() {
      //seqs = g_seqs;
      batch_num = g_batch_num;
      seqs = (bseq1_t*)malloc(batch_num*sizeof(bseq1_t));

      // duplicate bseq1_t so that tests can make modifications
      dup_bseq1_t(seqs, g_seqs, batch_num);
    }

    virtual void TearDown() {
      for (int i = 0; i < batch_num; i++) {
        if (seqs[i].name) free(seqs[i].name);
        if (seqs[i].comment) free(seqs[i].comment);
        if (seqs[i].seq) free(seqs[i].seq);
        if (seqs[i].qual) free(seqs[i].qual);
      }
      free(seqs);
    }

    bseq1_t* seqs;
    int batch_num;
};

class UtilTests : public BaseTests {
  ;
};

class PipelineTests : public BaseTests {
  ;
};

class BamTests : public ::testing::Test {
  ;
};

class BucketTest : public ::testing::Test {
  ;
};

class FPGATests;

#ifdef USE_MPI
class MPITests : public ::testing::Test {
 protected:
  virtual void SetUp() {
    if (MPI::COMM_WORLD.Get_size() > 1) {
      MPI::COMM_WORLD.Barrier();
    }
  }
  MPILink link_;
};

class ChannelTests : public MPITests {
  ;
};
#endif

inline std::string get_absolute_path(std::string path) {
  boost::filesystem::wpath file_path(path);
  if (boost::filesystem::exists(file_path)) {
    return boost::filesystem::canonical(file_path).string();
  }
  else {
    if (file_path.is_absolute()) {
      return path;
    }
    else {
      file_path = boost::filesystem::current_path() / file_path;
      return file_path.string();
    }
  }
}

inline std::string get_bin_dir() {
  std::stringstream ss;
  ss << "/proc/" << getpid() << "/exe";

  boost::filesystem::wpath bin_path(get_absolute_path(ss.str()));
  return bin_path.parent_path().string();
}

static inline void check_bam_core(bam1_core_t &base, bam1_core_t &test) {
    EXPECT_EQ(base.tid, test.tid);
    EXPECT_EQ(base.pos, test.pos);

    EXPECT_EQ(base.bin, test.bin);
    EXPECT_EQ(base.qual, test.qual);
    EXPECT_EQ(base.l_qname, test.l_qname);

    EXPECT_EQ(base.flag, test.flag);
    EXPECT_EQ(base.n_cigar, test.n_cigar);

    EXPECT_EQ(base.l_qseq, test.l_qseq);
    EXPECT_EQ(base.mtid, test.mtid);
    EXPECT_EQ(base.mpos, test.mpos);
    EXPECT_EQ(base.isize, test.isize);
}

static inline void check_bam(bam1_t& base, bam1_t& test) {
  check_bam_core(base.core, test.core);
  ASSERT_EQ(base.l_data, test.l_data);
  for (int i = 0; i < base.l_data; i++) {
    EXPECT_EQ(test.data[i], test.data[i]);
  }
}

static inline void check_bams(bams_t& base, bams_t& test) {
  ASSERT_EQ(base.l, test.l); 
  //EXPECT_EQ(1, base.l);
  for (int i = 0; i < base.l; i ++) {
    check_bam(*base.bams[i], *test.bams[i]);
  }
}

static inline void check_bseq(bseq1_t& base, bseq1_t& test) {
  ASSERT_EQ(base.l_seq, test.l_seq);
  // seq may not be encoded as string
  for (int i = 0; i < base.l_seq; i++) {
    EXPECT_EQ(base.seq[i], test.seq[i]);
  }
  
  EXPECT_STREQ(base.name, test.name);
  EXPECT_STREQ(base.comment, test.comment);
  EXPECT_STREQ(base.qual, test.qual);
}

static inline void check_mem_chain(mem_chain_t& c1, mem_chain_t& c2) {
  // check result match for mem_chain_t
  EXPECT_EQ(c1.n, c2.n);
  EXPECT_EQ(c1.m, c2.m);
  EXPECT_EQ(c1.first, c2.first);
  EXPECT_EQ(c1.rid, c2.rid);
  EXPECT_EQ(c1.w, c2.w);
  EXPECT_EQ(c1.kept, c2.kept);
  EXPECT_EQ(c1.is_alt, c2.is_alt);
  EXPECT_EQ(c1.frac_rep, c2.frac_rep);
  EXPECT_EQ(c1.pos, c2.pos);
}

static void check_mem_alnreg(mem_alnreg_t& a1, mem_alnreg_t& a2) {
  EXPECT_EQ(a1.rb, a2.rb);
  EXPECT_EQ(a1.re, a2.re);
  EXPECT_EQ(a1.qb, a2.qb);
  EXPECT_EQ(a1.qe, a2.qe);
	EXPECT_EQ(a1.rid, a2.rid);
	EXPECT_EQ(a1.score, a2.score);
	EXPECT_EQ(a1.truesc, a2.truesc);
	EXPECT_EQ(a1.sub, a2.sub);
	EXPECT_EQ(a1.alt_sc, a2.alt_sc);
	EXPECT_EQ(a1.csub, a2.csub);
	EXPECT_EQ(a1.sub_n, a2.sub_n);
	EXPECT_EQ(a1.w, a2.w);
	EXPECT_EQ(a1.seedcov, a2.seedcov);
	EXPECT_EQ(a1.secondary, a2.secondary);
	EXPECT_EQ(a1.secondary_all, a2.secondary_all);
	EXPECT_EQ(a1.seedlen0, a2.seedlen0);
  EXPECT_EQ(a1.n_comp, a2.n_comp);
  EXPECT_EQ(a1.is_alt, a2.is_alt);
	EXPECT_EQ(a1.frac_rep, a2.frac_rep);
	EXPECT_EQ(a1.hash, a2.hash);
}

#endif
