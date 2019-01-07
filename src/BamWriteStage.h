#include "Pipeline.h"

class BamWriteStage
: public kestrelFlow::MapStage<
      BamRecord, int, COMPUTE_DEPTH, COMPUTE_DEPTH>
{
 public:
  BamWriteStage(
      std::string bam_dir, 
      bam_hdr_t* h = NULL,
      int n = 1): 
    kestrelFlow::MapStage<
      BamRecord, int, COMPUTE_DEPTH, COMPUTE_DEPTH>(n), 
      bam_dir_(bam_dir), h_(h)
    {;} 

  int compute(BamRecord const & input);

 private:
  std::string  bam_dir_;
  bam_hdr_t*   h_; 
};
