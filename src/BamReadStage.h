#ifndef BAMREADSTAGE_H
#define BAMREADSTAGE_H

#include "Pipeline.h"

class BamReadStage
: public kestrelFlow::MapStage<
      int, BamRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>
{
 public:
  BamReadStage(
      std::string bam_dir, 
      bam_hdr_t* h = NULL,
      int n = 1): 
    kestrelFlow::MapStage<
      int, BamRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>(n), 
      bam_dir_(bam_dir), h_(h)
    {;} 

  BamRecord compute(int const & id);

 private:
  boost::mutex hdr_lock_;

  std::string  bam_dir_;
  bam_hdr_t*   h_;
};

#endif
