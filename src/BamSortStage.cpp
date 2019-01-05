#include "htslib/ksort.h"
#include "BamSortStage.h"

BamRecord BamSortStage::compute(BamRecord const & input) {
  BamRecord output = input;
  sort_bams(output.size, output.bams);
  return output;
}
