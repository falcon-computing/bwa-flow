#ifndef FPGA_TESTS_H
#define FPGA_TESTS_H

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
#include "TestCommon.h"

#include "FPGAAgent.h"
#include "FPGAPipeline.h"

namespace kestrelFlow {
class FPGATests : public BaseTests {
  public:
    void SetUp() {
      BaseTests::SetUp();
      bit_path_ = get_bin_dir() + "/../data/bit.awsxclbin";
      if (std::getenv("BWA_PAC_PATH")) {
        pac_path_ = std::getenv("BWA_PAC_PATH");
      }
      else {
        pac_path_ = "";
      }
    }

    std::string bit_path_;
    std::string pac_path_;
};
}

#endif
