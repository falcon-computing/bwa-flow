#include "Pipeline.h"
#include "samblaster.h"
#include <boost/thread/mutex.hpp>

// mapStage version
class MarkDupPartStage: public kestrelFlow::MapPartitionStage<
  BamsRecord, BamsRecord, INPUT_DEPTH, OUTPUT_DEPTH> {
public:
  MarkDupPartStage(ktp_aux_t* auxx = NULL):kestrelFlow::MapPartitionStage<
    BamsRecord, BamsRecord, INPUT_DEPTH, OUTPUT_DEPTH>(1, false){
      InitializeState(auxx);
      aux = auxx;
    }
  ~MarkDupPartStage(){
    // deleteState(state);
  } 
  void compute(int wid);
private:
  void InitializeState(ktp_aux_t* auxx); 
  state_t* state;
  ktp_aux_t* aux;
  boost::mutex mtx_;
};
