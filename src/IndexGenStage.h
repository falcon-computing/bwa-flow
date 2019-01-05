#include "Pipeline.h"

class IndexGenStage : public kestrelFlow::SourceStage<int, INPUT_DEPTH> {
 public:
  IndexGenStage(int n): 
    kestrelFlow::SourceStage<int, INPUT_DEPTH>(), num_ids_(n) {;} 
  void compute();
 private:
  int num_ids_;
};
