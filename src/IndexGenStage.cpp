#include "IndexGenStage.h"

void IndexGenStage::compute() {
  for (int i = 0; i < num_ids_; i++) {
    pushOutput(i);
  }
}
