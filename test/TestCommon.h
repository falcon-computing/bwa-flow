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

// Global variables
extern char *bwa_pg;
extern gzFile fp_idx, fp2_read2;
extern void *ko_read1, *ko_read2;
extern ktp_aux_t *aux;

class CommTests : public ::testing::Test {
  protected:
    virtual void SetUp() {
      // Read from file input, get mem_chains
      seqs = bseq_read(aux->actual_chunk_size, &batch_num, aux->ks, aux->ks2);
    }

    virtual void TearDown() {
      for (int i = 0; i < batch_num; i++) {
        free(seqs[i].name); 
        free(seqs[i].comment);
        free(seqs[i].seq); 
        free(seqs[i].qual); 
      }
      free(seqs);
    }

    bseq1_t* seqs;
    int batch_num;
};

#endif
