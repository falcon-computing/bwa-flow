#include "IndexGenStage.h"
#include "BamReadStage.h"
#include "BamSortStage.h"
#include "ReorderAndWriteStage.h"
#include "Pipeline.h"
#include "TestCommon.h"

namespace kf = kestrelFlow;

const int buckets = 4;
const int records_per_bucket = 16;
const std::string basedir = "./data/bucket-bams/";
const std::string baseline_file = "./data/bucket-bams/baseline.bam";
const std::string output_file = "./data/bucket-bams/output.bam";

TEST_F(BamTests, TestBamReadStage) {
  BamReadStage s(basedir);

  for (int i = 0; i < buckets; i++) {
    try {
      BamRecord b = s.compute(i);
      ASSERT_EQ(i, b.id);
      ASSERT_EQ(records_per_bucket, b.size);
    }
    catch (...) {
      FAIL() << "should not catch exceptions";
    }
  }

  try {
    BamRecord b = s.compute(buckets);
    FAIL() << "missing file should trigger an exception";
  }
  catch (...) {
    ;
  }
}

TEST_F(BamTests, TestBamSortStage) {
  BamReadStage s1(basedir);
  BamSortStage s2();

}

TEST_F(BamTests, TestReorderAndWriteStage) {
  samFile * fp = hts_open(baseline_file.c_str(), "r");
  bam_hdr_t * h = sam_hdr_read(fp);

  {
  BamReadStage s1(basedir);
  BamSortStage s2;
  ReorderAndWriteStage s3(output_file, h);

  kf::Pipeline p(3, 1);
  p.addStage(0, &s1);
  p.addStage(1, &s2);
  p.addStage(2, &s3);

  p.start();

  kf::Queue<int>* q = static_cast<kf::Queue<int>*>
                      (p.getInputQueue().get());
  q->push(3);
  q->push(0);
  q->push(1);
  q->push(2);

  p.finalize();
  p.wait();
  }

  // read bam files and make sure it's sorted
  std::stringstream ss;
  ss << "diff " << output_file << " " << baseline_file;
  int ret = system(ss.str().c_str());
  ASSERT_EQ(ret, 0);
  
}
