#include "Pipeline.h"

class ReorderAndWriteStage
: public kestrelFlow::SinkStage<BamRecord, COMPUTE_DEPTH>
{
  public:
    ReorderAndSortStage(std::string str, bam_hdr_t * h):kestrelFlow::SinkStage
                  <BamRecord, COMPUTE_DEPTH>(1, false), output_dir_(output_dir), h_(h) {
      fout_ = sam_open(str.c_str(), "wb");
      sam_hdr_write(fout_, h_);
    }
    ~ReorderAndWriteStage() {
      sam_close(fout_);
    }
    void compute();
  private:
    std::string output_dir_;
    samFile * fout_;
    bam_hdr_t * h_;
};