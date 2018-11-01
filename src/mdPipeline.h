#include "Pipeline.h"
#include "samblaster.h"
#include "config.h"

#include "boost/thread/mutex.hpp"

class Markdup: public kestrelFlow::MapStage<
  SeqsRecord, SeqsRecord, INPUT_DEPTH, OUTPUT_DEPTH> {
public:
  Markdup(int n=1, ktp_aux_t* auxx = NULL):kestrelFlow::MapStage<
    SeqsRecord, SeqsRecord, INPUT_DEPTH, OUTPUT_DEPTH>(n, false){
      InitializeState(auxx);
      aux = auxx;
    }

  SeqsRecord compute(SeqsRecord const & record);
private:
  void InitializeState(ktp_aux_t* auxx); 
  state_t* state;
  ktp_aux_t* aux;

  boost::mutex mtx_;
};
