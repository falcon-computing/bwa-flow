#include "BamSortStage.h"
#include "config.h"
#include "htslib/ksort.h"
#include "util.h"

BamRecord BamSortStage::compute(BamRecord const & input) {
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started BamSort()";

  BamRecord output = input;

  if (!FLAGS_disable_sort) {
    uint64_t start_ts = getUs();
    sort_bams(output.size, output.bams);

    DLOG_IF(INFO, VLOG_IS_ON(1)) << "sorting " << output.size
      << " records took " << getUs() - start_ts << " us";
  }

  // write serialized BAM to a buffer first
  // NOTE: not sure how much benefits so far
  uint64_t start_ts = getUs(); 
  output.fbuf = new BamFileBuffer(0);
  for (int i = 0; i < output.size; i++) {
    if (output.fbuf->write(output.bams[i]) < 0) {
      throw std::runtime_error("cannot convert bam file");
    }
    bam_destroy1(output.bams[i]);
  }
  free(output.bams);

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "serializing " << output.size
    << " records took " << getUs() - start_ts << " us";

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished Bamsort()";
  return output;
}
