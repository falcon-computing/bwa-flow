#include "Pipeline.h"
#include "samblaster.h"
#include <boost/thread/mutex.hpp>

// mapStage version
class MarkDupStage: public kestrelFlow::MapStage<
  SeqsRecord, SeqsRecord, INPUT_DEPTH, OUTPUT_DEPTH> {
public:
  MarkDupStage(int n = 1, ktp_aux_t* auxx = NULL):kestrelFlow::MapStage<
    SeqsRecord, SeqsRecord, INPUT_DEPTH, OUTPUT_DEPTH>(n){
      InitializeState(auxx);
      aux = auxx;
    }
  ~MarkDupStage() {
    // deleteState(state);
  } 
  SeqsRecord compute(SeqsRecord const & input);
private:
  void InitializeState(ktp_aux_t* auxx); 
  state_t* state;
  ktp_aux_t* aux;
  boost::mutex mtx_;
};
