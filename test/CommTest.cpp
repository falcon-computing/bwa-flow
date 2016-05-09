#include "bwa_wrapper.h"
#include "Pipeline.h"
#include "TestCommon.h"

TEST_F(CommTests, InputSerializationTest) {

  // Read from file input, get mem_chains
  int batch_num = 0;
  bseq1_t *seqs = bseq_read(aux.actual_chunk_size, 
      &batch_num, aux.ks, aux.ks2);

  SeqsRecord record;
  record.start_idx = 0x1027;    // use a magic number for testing
  record.batch_num = batch_num;
  record.seqs = seqs;

  SeqsDispatch send_stage;
  SeqsReceive  recv_stage;

  std::string ser_data = send_stage.serialize(&record);

  ASSERT_GT(ser_data.size(), 0);

  SeqsRecord record_test = recv_stage.deserialize(
      ser_data.c_str(), ser_data.length());

  ASSERT_EQ(record.start_idx, record_test.start_idx);
  ASSERT_EQ(record.batch_num, record_test.batch_num);

  for (int i = 0; i < batch_num; i++) {
    bseq1_t* seq_base = &seqs[i];
    bseq1_t* seq_test = &record_test.seqs[i];

    ASSERT_EQ(seq_base->l_seq, seq_test->l_seq);
    ASSERT_STREQ(seq_base->name, seq_test->name);
    ASSERT_STREQ(seq_base->comment, seq_test->comment);
    ASSERT_STREQ(seq_base->seq, seq_test->seq);
    ASSERT_STREQ(seq_base->qual, seq_test->qual);
  }

  freeSeqs(seqs, batch_num);
  freeSeqs(record_test.seqs, batch_num);
}
