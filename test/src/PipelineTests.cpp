#include "kflow/Pipeline.h"
#include "Pipeline.h"
#include "TestCommon.h"

TEST_F(PipelineTests, NewSeqRead) {
  // KseqsRead       kread_stage;
  // KseqsToBseqs    k2b_stage;

  // kestrelFlow::Pipeline pipeline(2, 1);
  // pipeline.addStage(0, &kread_stage);
  // pipeline.addStage(1, &k2b_stage);
  
  // TODO: need to re-initialize aux for file read
}

TEST_F(PipelineTests, Seq2BamsCompute) {
  int test_num = batch_num > 1024 ? 1024 : batch_num;

  SeqsRecord input;
  input.start_idx = 0;
  input.batch_num = test_num;
  input.seqs = seqs;

  // compute original results
  bseq1_t* base = bwa_mem(seqs, test_num);
  
  // test pipeline with single thread
  SeqsToChains    stage_1(1);
  ChainsToRegions stage_2(1);
  RegionsToSam    stage_3(1);

  ChainsRecord  inter_1 = stage_1.compute(input);
  RegionsRecord inter_2 = stage_2.compute(inter_1);
  SeqsRecord    output  = stage_3.compute(inter_2);

  for (int i = 0; i < test_num; i ++) {
    check_bams(base[i].bams[0], seqs[i].bams[0]);

    // free bams_t
    for (int j = 0; j < base[i].bams->l; j++) {
      bam_destroy1(seqs[i].bams->bams[j]);
      bam_destroy1(base[i].bams->bams[j]);
    }
    // avoid double-free
    seqs[i].name = 0;
    seqs[i].comment = 0;
    seqs[i].seq = 0;
    seqs[i].qual = 0;

    free(base[i].bams);
    free(seqs[i].bams);
  }

  free(base);
}

