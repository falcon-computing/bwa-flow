#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include <unistd.h>
#include <dlfcn.h>

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

#ifdef USE_MPI
class MPITests : public ::testing::Test {
 protected:
  MPILink link_;
};

class ChannelTests : public MPITests {
  ;
};
#endif

#endif
