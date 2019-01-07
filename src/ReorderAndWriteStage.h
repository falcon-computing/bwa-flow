#ifndef REORDERANDWRITESTAGE_H
#define REORDERANDWRITESTAGE_H

#include "Pipeline.h"

class ReorderAndWriteStage
: public kestrelFlow::SinkStage<BamRecord, COMPUTE_DEPTH>
{
 public:
  ReorderAndWriteStage(std::string output, bam_hdr_t * h):
    kestrelFlow::SinkStage<BamRecord, COMPUTE_DEPTH>(1, false), 
    output_(output), h_(h) 
  {
    fout_ = hts_open(output_.c_str(), "wb0");
    sam_hdr_write(fout_, h_);
  }

  ~ReorderAndWriteStage() {
    sam_close(fout_);
  }
  void compute(int wid);

 private:
  std::string output_;
  samFile * fout_;
  bam_hdr_t * h_;
};

#endif
