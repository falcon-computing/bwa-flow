#include "Pipeline.h"
#include "samblaster.h"
#include <boost/thread/mutex.hpp>

// mapStage version
class MarkDup: public kestrelFlow::MapStage<
  SeqsRecord, SeqsRecord, INPUT_DEPTH, OUTPUT_DEPTH> {
public:
  MarkDup(int n = 1, ktp_aux_t* auxx = NULL):kestrelFlow::MapStage<
    SeqsRecord, SeqsRecord, INPUT_DEPTH, OUTPUT_DEPTH>(n){
      InitializeState(auxx);
      aux = auxx;
    }
  ~MarkDup() {
    // deleteState(state);
  } 
  SeqsRecord compute(SeqsRecord const & input);
private:
  void InitializeState(ktp_aux_t* auxx); 
  state_t* state;
  ktp_aux_t* aux;
  boost::mutex mtx_;
};
