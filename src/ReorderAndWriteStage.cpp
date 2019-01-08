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

      bgzf_write(fout_->fp.bgzf, 
          record.fbuf->get_data(), 
          record.fbuf->get_size());

      delete record.fbuf;
      
      DLOG(INFO) << "Wrtting " << record.fbuf->get_size() / 1024 << " kB took " 
        << getUs() - start_ts << " us";

    }
  }
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished ReorderAndWrite()";
}
