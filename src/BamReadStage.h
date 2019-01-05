#include "Pipeline.h"

class BamReadStage
: public kestrelFlow::MapStage<
      int, BamRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>
{
 public:
  BamReadStage(std::string temp_dir, int n=1): kestrelFlow::MapStage<
      int, BamRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>(n), temp_dir_(temp_dir) {;} 

  BamRecord compute(int const & id);
 private:
  std::string temp_dir_;
};
