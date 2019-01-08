#include "IndexGenStage.h"

void IndexGenStage::compute() {
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started IndexGenStage()";
  for (int i = 0; i < num_ids_; i++) {
    pushOutput(i);
  }
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished IndexGenStage()";
}
