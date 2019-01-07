#include "Pipeline.h"

class BamWriteStage
: public kestrelFlow::MapStage<
      BamRecord, int, COMPUTE_DEPTH, COMPUTE_DEPTH>
{
 public:
  BamWriteStage(
      int num_parts,
      std::string bam_dir, 
      std::string output_path, 
      bam_hdr_t* h = NULL,
      int n = 1);

  ~BamWriteStage();

  int compute(BamRecord const & input);

 private:
  int         num_parts_;
  std::string bam_dir_;
  std::string output_path_;
  bam_hdr_t*  h_; 
};
