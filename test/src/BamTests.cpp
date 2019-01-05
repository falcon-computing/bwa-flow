#include "IndexGenStage.h"
#include "BamReadStage.h"
#include "ReorderAndWriteStage.h"
#include "Pipeline.h"
#include "TestCommon.h"

const int buckets = 4;
const int records_per_bucket = 16;
const std::string basedir = "./bucket-bams/";
const std::string baseline_file = "./bucket-bams/baseline.bam";
const std::string output_file = "./bucket-bams/output.bam";

TEST_F(BamTests, TestBamReadStage) {
  BamReadStage s(basedir);

  for (int i = 0; i < buckets; i++) {
    try {
      BamRecord b = s.compute(i);
      ASSERT_EQ(records_per_bucket, b.size);
    }
    catch (..) {
      FAIL() << "should not catch exceptions";
    }
  }

  try {
    BamRecord b = s.compute(buckets);
    FAIL() << "missing file should trigger an exception";
  }
  catch (..) {
    ;
  }
}

TEST_F(BamTests, TestBamSortStage) {
  BamReadStage s1(basedir);
  BamSortStage s2;

}

TEST_F(BamTests, TestReorderAndWriteStage) {
  BamReadStage s1(basedir);
  BamSortStage s2;
  std::vector<BamRecord> records(buckets);

  for (int i = 0; i < buckets; i++) {
    try {
      records[i] = s2.compute(s1.compute(i));
    }
    catch (..) {
      FAIL() << "should not catch exceptions";
    }
  }

  /*
  {
  ReorderAndWriteStage s3();
  s3.push(records[3]);
  s3.push(records[0]);
  s3.push(records[1]);
  s3.push(records[2]);
  s3.finalize();
  }

  // read bam files and make sure it's sorted
  samFile* test = hts_open(output_file
  
  */
}
