#include "htslib/ksort.h"
#include "BamSortStage.h"

BamRecord BamSortStage::compute(BamRecord const & input) {
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started BamSort()";
  BamRecord output = input;
  sort_bams(output.size, output.bams);
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished Bamsort()";
  //if (input.size > 0) {
  //DLOG(INFO) <<"???" << output.bams[0]->core.tid;
  //}
  return output;
}
