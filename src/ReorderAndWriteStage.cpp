#include "ReorderAndWriteStage.h"

void ReorderAndWriteStage::compute(){
  uint32_t n_processed = 0;
  std::unordered_map<uint32_t, bam1_t**> record_buf;
  std::vector<SeqsRecord> output_records;
  int max_batch_records = FLAGS_max_batch_records;
  int batch_records = 0;
  int bam_buffer_order = 0;
  try {
   while (true) {
     BamRecord input;
     SeqsRecord record;
     bool ready = this->getInput(input);
     while (!this->isFinal() && !ready) {
       boost::this_thread::sleep_for(boost::chrono::microseconds(10));
       ready = this->getInput(input);
     }
     if (!ready) {
       break;
     }
     record_buf[input.id] = input;
     while (record_buf.count(n_processed)) {
         record = record_buf[n_processed];
         record_buf.erase(n_processed);
         n_processed += 1;

         for (int i = 0; i < input.size; i++) {
           sam_wirte1(fout_, h_, input.bams[i]);
           bam_destroy1(input.bams[i]);
         }

         free(input.bams);
     }
   }
  }
}
