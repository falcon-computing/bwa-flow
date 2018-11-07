#include "Pipeline.h"
#include "samblaster.h"

class MarkDup: public kestrelFlow::MapPartitionStage<
  SeqsRecord, SeqsRecord, INPUT_DEPTH, OUTPUT_DEPTH> {
public:
  MarkDup(ktp_aux_t* auxx = NULL):kestrelFlow::MapPartitionStage<
    SeqsRecord, SeqsRecord, INPUT_DEPTH, OUTPUT_DEPTH>(1, false){
      InitializeState(auxx);
      aux = auxx;
    }
  ~MarkDup(){
    // deleteState(state);
  } 
  void compute(int wid);
private:
  void InitializeState(ktp_aux_t* auxx); 
  state_t* state;
  ktp_aux_t* aux;
};
