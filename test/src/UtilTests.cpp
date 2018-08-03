#include "bwa_wrapper.h"
#include "Pipeline.h"
#include "util.h"
#include "TestCommon.h"

TEST_F(UtilTests, SeqTest) {
  int test_num = batch_num > 4096 ? 4096 : batch_num;
  uint64_t start_ts = getUs();
  std::stringstream ss;
  for (int i = 0; i < test_num; i++) {
    serialize(ss, seqs[i]);
  }
  ASSERT_LT(0, ss.str().size());

  VLOG(1) << test_num << " seqs (for chain record) has "
          << ss.str().size() << " bytes, "
          << "serialization takes "
          << getUs() - start_ts << " us";
  
  start_ts = getUs();
  bseq1_t* seqs_test = (bseq1_t*)malloc(test_num*sizeof(bseq1_t));
  for (int i = 0; i < test_num; i++) {
    deserialize(ss, seqs_test[i]);
  }
  VLOG(1) << "deserialization takes "
          << getUs() - start_ts << " us";

  // skip testing
}

TEST_F(UtilTests, MemChainTest) {
  int test_num = batch_num > 32 ? 32 : batch_num;
  // here only test the first 32 to save time
  for (int i = 0; i < test_num; i++) {
    mem_chain_v chains = mem_seq2chain_wrapper(aux, seqs);

    std::stringstream ss;
    serialize(ss, chains);

    ASSERT_LT(0, ss.str().size());

    mem_chain_v chains_test;
    deserialize(ss, chains_test);

    ASSERT_EQ(chains.n, chains_test.n);
    ASSERT_EQ(chains.m, chains_test.m);

    for (int j = 0; j < chains.n; j++) {
      check_mem_chain(chains.a[j], chains_test.a[j]);
    }
    
    // free chains
    for (int j = 0; j < chains.n; j++) {
      free(chains.a[j].seeds);
      free(chains_test.a[j].seeds);
    }
    free(chains.a);
    free(chains_test.a);
  }
}

TEST_F(UtilTests, MemAlnregTest) {
  int test_num = batch_num > 32 ? 32 : batch_num;
  for (int i = 0; i < test_num; i++) {
    mem_chain_v chains = mem_seq2chain_wrapper(aux, seqs);
    mem_alnreg_v alnregs;
    kv_init(alnregs);
    for (int j = 0; j < chains.n; j++) {
      mem_chain2aln(
          aux->opt, 
          aux->idx->bns, 
          aux->idx->pac,
          seqs[i].l_seq,
          (uint8_t*)seqs[i].seq,
          &chains.a[j],
          &alnregs);
    }

    // serialize mem_alnreg_v
    std::stringstream ss;
    serialize(ss, alnregs);

    ASSERT_LT(0, ss.str().size());

    mem_alnreg_v alnregs_test;
    deserialize(ss, alnregs_test);

    ASSERT_EQ(alnregs.n, alnregs_test.n);
    ASSERT_EQ(alnregs.m, alnregs_test.m);

    for (int j = 0; j < alnregs.n; j++) {
      check_mem_alnreg(alnregs.a[j], alnregs_test.a[j]);
    }
  }
}

TEST_F(UtilTests, MemChainVFullTest) {
  int test_num = batch_num > 4096 ? 4096 : batch_num;

  mem_chain_v* chains_base = (mem_chain_v*)malloc(sizeof(mem_chain_v)*test_num);
  mem_chain_v* chains_test = (mem_chain_v*)malloc(sizeof(mem_chain_v)*test_num);
  for (int i = 0; i < test_num; i++) {
    chains_base[i] = mem_seq2chain_wrapper(aux, seqs);
  }

  // start serialization
  std::stringstream ss;
  uint64_t start_ts = getUs(); 
  for (int i = 0; i < test_num; i++) {
    serialize(ss, chains_base[i]);
  }
  VLOG(1) << test_num << " chains has "
          << ss.str().size() << " bytes, "
          << "serialization takes "
          << getUs() - start_ts << " us";

  // start deserialization
  start_ts = getUs();
  for (int i = 0; i < test_num; i++) {
    deserialize(ss, chains_test[i]);
  }
  VLOG(1) << "deserialization takes "
          << getUs() - start_ts << " us";

  // verify
  for (int i = 0; i < test_num; i++) {
    ASSERT_EQ(chains_base[i].n, chains_test[i].n);
    ASSERT_EQ(chains_base[i].m, chains_test[i].m);

    for (int j = 0; j < chains_base[i].n; j++) {
      check_mem_chain(chains_base[i].a[j], chains_test[i].a[j]);
      // free chains
      free(chains_base[i].a[j].seeds);
      free(chains_test[i].a[j].seeds);
    }
    free(chains_base[i].a);
    free(chains_test[i].a);
  }
  free(chains_base);
  free(chains_test);
}

TEST_F(UtilTests, MemAlnregVFullTest) {
  int test_num = batch_num > 4096 ? 4096 : batch_num;

  mem_alnreg_v* aln_base = (mem_alnreg_v*)malloc(test_num*sizeof(mem_alnreg_v));
  mem_alnreg_v* aln_test = (mem_alnreg_v*)malloc(test_num*sizeof(mem_alnreg_v)); 

  for (int i = 0; i < test_num; i++) {
    mem_chain_v chains = mem_seq2chain_wrapper(aux, seqs);
    kv_init(aln_base[i]);
    for (int j = 0; j < chains.n; j++) {
      mem_chain2aln(
          aux->opt, 
          aux->idx->bns, 
          aux->idx->pac,
          seqs[i].l_seq,
          (uint8_t*)seqs[i].seq,
          &chains.a[j],
          &aln_base[i]);
    }
  }

  // start serialization
  std::stringstream ss;
  uint64_t start_ts = getUs(); 
  for (int i = 0; i < test_num; i++) {
    serialize(ss, aln_base[i]);
  }
  VLOG(1) << test_num << " regions has "
          << ss.str().size() << " bytes, "
          << "serialization takes "
          << getUs() - start_ts << " us";

  // start deserialization
  start_ts = getUs();
  for (int i = 0; i < test_num; i++) {
    deserialize(ss, aln_test[i]);
  }
  VLOG(1) << "deserialization takes "
          << getUs() - start_ts << " us";

  // verify
  for (int i = 0; i < test_num; i++) {
    ASSERT_EQ(aln_base[i].n, aln_test[i].n);
    ASSERT_EQ(aln_base[i].m, aln_test[i].m);

    for (int j = 0; j < aln_base[i].n; j++) {
      check_mem_alnreg(aln_base[i].a[j], aln_test[i].a[j]);
    }
    free(aln_base[i].a);
    free(aln_test[i].a);
  }
  free(aln_base);
  free(aln_test);
}

TEST_F(UtilTests, SeqsRecordTest) {
  int test_num = batch_num > 4096 ? 4096 : batch_num;

  SeqsRecord record;
  record.start_idx = 0x1027;    // use a magic number for testing
  record.batch_num = test_num;
  record.seqs = seqs;

  uint64_t start_ts = getUs();
  std::string ser_data = serialize(record);

  ASSERT_GT(ser_data.size(), 0);

  VLOG(1) << test_num << " seqs has "
          << ser_data.size() << " bytes, "
          << "serialization takes "
          << getUs() - start_ts << " us";

  start_ts = getUs();
  SeqsRecord record_test;
  deserialize(ser_data.c_str(), ser_data.length(), record_test);

  VLOG(1) << "deserialization takes "
          << getUs() - start_ts << " us";

  ASSERT_EQ(record.start_idx, record_test.start_idx);
  ASSERT_EQ(record.batch_num, record_test.batch_num);

  for (int i = 0; i < test_num; i++) {
    check_bseq(seqs[i], record_test.seqs[i]); 
  }

  freeSeqs(record_test.seqs, test_num);
}

TEST_F(UtilTests, ChainsRecordTest) {
  int test_num = batch_num > 4096 ? 4096 : batch_num;

  // compute mem_chain_v
  mem_chain_v* chains = (mem_chain_v*)malloc(sizeof(mem_chain_v)*test_num);
  for (int i = 0; i < test_num; i++) {
    chains[i] = mem_seq2chain_wrapper(aux, seqs);
  }

  // construct record
  ChainsRecord record;
  record.start_idx = 0x1027;    // use a magic number for testing
  record.batch_num = test_num;
  record.seqs = seqs;
  record.chains = chains;

  uint64_t start_ts = getUs();

  // serialization
  std::string ser_data = serialize(record);
  ASSERT_GT(ser_data.size(), 0);

  VLOG(1) << test_num << " chains has "
          << ser_data.size() << " bytes, "
          << "serialization takes "
          << getUs() - start_ts << " us";

  start_ts = getUs();

  // deserialization
  ChainsRecord test;
  deserialize(ser_data.c_str(), ser_data.length(), test);

  VLOG(1) << "deserialization takes "
          << getUs() - start_ts << " us";

  // check results
  ASSERT_EQ(record.start_idx, test.start_idx);
  ASSERT_EQ(record.batch_num, test.batch_num);

  for (int i = 0; i < test_num; i++) {
    check_bseq(seqs[i], test.seqs[i]); 

    ASSERT_EQ(record.chains[i].n, test.chains[i].n);
    for (int j = 0; j < record.chains[i].n; j++) {
      check_mem_chain(record.chains[i].a[j], test.chains[i].a[j]); 
    }
  }

  // clean up
  freeSeqs(test.seqs, test_num);
  freeChains(test.chains, test_num);
}

TEST_F(UtilTests, RegionsRecordTest) {
  int test_num = batch_num > 4096 ? 4096 : batch_num;

  // compute mem_alnreg_v
  mem_alnreg_v* alnreg = (mem_alnreg_v*)malloc(test_num*sizeof(mem_alnreg_v));
  for (int i = 0; i < test_num; i++) {
    mem_chain_v chains = mem_seq2chain_wrapper(aux, seqs);
    kv_init(alnreg[i]);
    for (int j = 0; j < chains.n; j++) {
      mem_chain2aln(
          aux->opt, 
          aux->idx->bns, 
          aux->idx->pac,
          seqs[i].l_seq,
          (uint8_t*)seqs[i].seq,
          &chains.a[j],
          alnreg+i);
    }
  }

  // construct record
  RegionsRecord record;
  record.start_idx = 0x1027;    // use a magic number for testing
  record.batch_num = test_num;
  record.seqs = seqs;
  record.alnreg = alnreg;

  uint64_t start_ts = getUs();

  // serialization
  std::string ser_data = serialize(record);
  ASSERT_GT(ser_data.size(), 0);

  VLOG(1) << test_num << " chains has "
          << ser_data.size() << " bytes, "
          << "serialization takes "
          << getUs() - start_ts << " us";

  start_ts = getUs();

  // deserialization
  RegionsRecord test;
  deserialize(ser_data.c_str(), ser_data.length(), test);

  VLOG(1) << "deserialization takes "
          << getUs() - start_ts << " us";

  // check results
  ASSERT_EQ(record.start_idx, test.start_idx);
  ASSERT_EQ(record.batch_num, test.batch_num);

  for (int i = 0; i < test_num; i++) {
    check_bseq(seqs[i], test.seqs[i]); 

    ASSERT_EQ(record.alnreg[i].n, test.alnreg[i].n);
    for (int j = 0; j < record.alnreg[i].n; j++) {
      check_mem_alnreg(record.alnreg[i].a[j], test.alnreg[i].a[j]); 
    }
  }

  freeSeqs(test.seqs, test_num);
  freeAligns(test.alnreg, test_num);
}
