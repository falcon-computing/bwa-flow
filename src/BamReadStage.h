#include <Pipeline.h>

class BamReadStage
: public kestrelFlow::MapStage<
      int, BamRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>
{
 public:
  BamReadStage(int n=1): kestrelFlow::MapStage<
      int, BamRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>(n) {;} 

  BamRecord compute(int id);
};
