#define TEST_FRIENDS_LIST \
          FRIEND_TEST(FPGATests, FPGATest); 

#include "TestCommon.h"
#include "FPGATests.h"

#include "bwa_wrapper.h"
#include "kflow/Stage.h"
#include "kflow/Queue.h"
#include "kflow/Pipeline.h"
#include "Pipeline.h"

#include "FPGAAgent.h"
#include "FPGAPipeline.h"

BWAOCLEnv* opencl_env;

namespace kestrelFlow {
TEST_F(FPGATests, PACTest) {
  if (pac_path_.empty()) {
    std::cout << "[   SKIP   ] "
              << "pac_path is missing"
              << std::endl;
    return;
  }

  char* pac;
  int64_t pac_size = BWAOCLEnv::get_full_pac(pac);
  int64_t pac_size_f;

  // load PAC file from fpga.pac
  FILE* pac_fpga = fopen(pac_path_.c_str(), "rb");
  fseek(pac_fpga, 0, SEEK_END);
  pac_size_f = ftell(pac_fpga);
  fseek(pac_fpga, 0, SEEK_SET);

  ASSERT_EQ(pac_size_f, pac_size+1);

  char* pac_2 = (char*)calloc(pac_size_f, sizeof(char));
  fread(pac_2, 1, pac_size_f-1, pac_fpga);
  fclose(pac_fpga);

  for (int64_t k = 0; k < pac_size_f-1; k++) {
    ASSERT_EQ(pac_2[k], pac[k]) << "k = " << k;
  }

  free(pac);
  free(pac_2);
}

TEST_F(FPGATests, FPGATest) {
  int test_num = batch_num > 1024 ? 1024 : batch_num;
  //int test_num = batch_num;

  try {
    FLAGS_fpga_path = bit_path_;
    opencl_env = new BWAOCLEnv();

    SeqsRecord input;
    input.start_idx = 0;
    input.batch_num = test_num;
    input.seqs = seqs;

    // compute base results
    bseq1_t* base = bwa_mem(seqs, test_num);

    // compute standard input/output using CPU
    Pipeline p(3, 1);
    SeqsToChains        stage_1(1);
    ChainsToRegionsFPGA stage_2(1);
    RegionsToSam        stage_3(1);

    p.addStage(0, &stage_1);
    p.addStage(1, &stage_2);
    p.addStage(2, &stage_3);

    // compute fpga results with 1 cpu and 1 fpga
    MegaPipe mp(1, 1);
    mp.addPipeline(&p);

    // input and output queue
    typedef Queue<SeqsRecord, INPUT_DEPTH> QIN;
    QIN* iq = static_cast<QIN*>(p.getInputQueue().get());
    if (!iq) { 
      FAIL() << "cannot case input queue";
    }

    // start FPGA computation
    DLOG(INFO) << "push one input to FPGA";
    iq->push(input);

    // after pipeline finishes seqs will contain the results
    mp.start();
    mp.finalize();
    mp.wait();

    // check results
    // the fpga results of regions record is unfiltered,
    // so need to compare the sam record instead
    for (int i = 0; i < test_num; i++) {
      check_bams(base[i].bams[0], seqs[i].bams[0]);

      // free bams_t
      for (int j = 0; j < base[i].bams->l; j++) {
        bam_destroy1(base[i].bams->bams[j]);
        bam_destroy1(seqs[i].bams->bams[j]);
      }
      // avoid double-free
      seqs[i].name = 0;
      seqs[i].comment = 0;
      seqs[i].seq = 0;
      seqs[i].qual = 0;

      free(base[i].bams);
      free(seqs[i].bams);

    }
    delete opencl_env;
  }
  catch (std::runtime_error &e) {
    std::cout << "[   SKIP   ] "
              << "cannot init FPGA because: " << e.what()
              << std::endl;
  }
}
}
