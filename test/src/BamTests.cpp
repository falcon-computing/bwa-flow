#include "IndexGenStage.h"
#include "BamReadStage.h"
#include "BamSortStage.h"
#include "BamFileBuffer.h"
#include "BamWriteStage.h"
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
  BamReadStage s1(basedir, aux->h);
  BamSortStage s2();

}

TEST_F(BamTests, TestBamWriteStage) {
  samFile * fp = hts_open(baseline_file.c_str(), "r");
  bam_hdr_t * h = sam_hdr_read(fp);

  {
  IndexGenStage s0(4);
  BamReadStage  s1(basedir, h);
  BamSortStage  s2;
  BamWriteStage s3(4, basedir, output_file, h);

  kf::Pipeline p(4, 1);
  p.addStage(0, &s0);
  p.addStage(1, &s1);
  p.addStage(2, &s2);
  p.addStage(3, &s3);

  kf::MegaPipe pipe(1, 0);
  pipe.addPipeline(&p);

  pipe.start();
  pipe.wait();
  }

  // read bam files and make sure it's sorted
  samFile * fp1 = hts_open(output_file.c_str(), "r");
  bam_hdr_t * h1 = sam_hdr_read(fp1);

  while (true) {
    bam1_t* b1 = bam_init1();
    bam1_t* b2 = bam_init1();
    int r1 = sam_read1(fp, h, b1);
    int r2 = sam_read1(fp1, h1, b2);
    if (r1 < 0 && r2 < 0) {
      break;
    }
    else if (r1 >= 0 && r2 >= 0) {
      check_bam(*b1, *b2);
    }
    else {
      FAIL() << "unmatched number of records";
    }
    bam_destroy1(b1);
    bam_destroy1(b2);
  }

  sam_close(fp);
  sam_close(fp1);
  
}

TEST_F(BamTests, TestReorderAndWriteStage) {
  samFile * fp = hts_open(baseline_file.c_str(), "r");
  bam_hdr_t * h = sam_hdr_read(fp);

  {
  IndexGenStage s0(4);
  BamReadStage s1(basedir, h);
  BamSortStage s2;
  ReorderAndWriteStage s3(output_file, h);

  kf::Pipeline p(4, 1);
  p.addStage(0, &s0);
  p.addStage(1, &s1);
  p.addStage(2, &s2);
  p.addStage(3, &s3);

  kf::MegaPipe pipe(1, 0);
  pipe.addPipeline(&p);

  pipe.start();
  pipe.wait();
  }

  // read bam files and make sure it's sorted
  samFile * fp1 = hts_open(output_file.c_str(), "r");
  bam_hdr_t * h1 = sam_hdr_read(fp1);

  while (true) {
    bam1_t* b1 = bam_init1();
    bam1_t* b2 = bam_init1();
    int r1 = sam_read1(fp, h, b1);
    int r2 = sam_read1(fp1, h1, b2);
    if (r1 < 0 && r2 < 0) {
      break;
    }
    else if (r1 >= 0 && r2 >= 0) {
      check_bam(*b1, *b2);
    }
    else {
      FAIL() << "unmatched number of records";
    }
    bam_destroy1(b1);
    bam_destroy1(b2);
  }

  sam_close(fp);
  sam_close(fp1);
  
}

TEST_F(BamTests, TestBamFileBuffer) {
  samFile* fp  = hts_open(baseline_file.c_str(), "r");
  bam_hdr_t* h = sam_hdr_read(fp);

  std::vector<bam1_t*> records(records_per_bucket*buckets);

  BamFileBuffer fb(fp->fp.bgzf->is_be);

  int i = 0;
  bam1_t* b = bam_init1();
  while (sam_read1(fp, h, b) > 0) {
    records[i++] = b;

    ASSERT_GT(fb.write(b), 0);
  }
  sam_close(fp);

  // write file to disk
  samFile* fout  = hts_open(output_file.c_str(), "wb");
  sam_hdr_write(fout, h);
  bgzf_write(fout->fp.bgzf, fb.get_data(), fb.get_size());
  sam_close(fout);

  // read file again
  samFile* fin  = hts_open(output_file.c_str(), "r");
  i = 0;
  while (sam_read1(fin, h, b) > 0) {
    check_bam(*b, *records[i]);
    bam_destroy1(records[i]);
    bam_destroy1(b);
    i++;
  }
  sam_close(fin);
}
