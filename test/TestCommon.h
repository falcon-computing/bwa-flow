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

class CommTests : public ::testing::Test {
  protected:
    CommTests() { }
};

// Global variables
extern char *bwa_pg;
extern gzFile fp_idx, fp2_read2;
extern void *ko_read1, *ko_read2;
extern ktp_aux_t aux;

#endif
