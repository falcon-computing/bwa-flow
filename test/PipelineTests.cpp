#include "bwa_wrapper.h"
#include "kflow/Pipeline.h"
#include "Pipeline.h"
#include "test/TestCommon.h"

bseq1_t* bwa_mem(bseq1_t* input, int batch_num) {

  // output bseq containing bams
  bseq1_t* seqs = (bseq1_t*)malloc(batch_num*sizeof(bseq1_t));
  dup_bseq1_t(seqs, input, batch_num);

  // temporary regions
  mem_alnreg_v* alnreg = new mem_alnreg_v[batch_num];

  for (int i = 0; i < batch_num; i++) {
    mem_chain_v chains = seq2chain(aux, &seqs[i]);
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
      mem_sam_pe(
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
      mem_reg2sam(
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
  EXPECT_EQ(1, base.l);
  for (int i = 0; i < base.l; i ++) {
    check_bam(*base.bams[i], *test.bams[i]);
  }
}

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

