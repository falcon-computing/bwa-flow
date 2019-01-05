#include "BamSortStage.h"

BamRecord BamSortStage::compute(BamRecord input) {
  ks_mergesort(sort, input.size, input.bams, 0);
  return input;
}
