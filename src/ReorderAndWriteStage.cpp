#include "ReorderAndWriteStage.h"
#include "util.h"

void ReorderAndWriteStage::compute(int wid) {
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started ReorderAndWrite()";
  uint32_t n_processed = 0;
  std::unordered_map<uint32_t, BamRecord> record_buf;

  while (true) {
    BamRecord input;
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
      BamRecord record = record_buf[n_processed];
      //record_buf.erase(n_processed);
      n_processed += 1;

      uint64_t start_ts = getUs();
      for (int i = 0; i < record.size; i++) {
        if (sam_write1(fout_, h_, record.bams[i]) < 0) {
          throw std::runtime_error("sam_write1() error");
        } 
        bam_destroy1(record.bams[i]);
      }
      DLOG(INFO) << "Wrtting " << record.size << " records took " 
        << getUs() - start_ts << " us";

      free(record.bams);
    }
  }
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished ReorderAndWrite()";
}
