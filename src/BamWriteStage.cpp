#include <sstream>
#include <iostream>
#include <iomanip>
#include "BamWriteStage.h"

int BamWriteStage::compute(BamRecord const & input) {
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started BamWrite()";
  
  int id = input.id;
  std::stringstream ss;
  ss << bam_dir_ << "/part-" << std::setw(6) << std::setfill('0')
     << id << "_sort.bam";
  samFile * fout = hts_open(ss.str().c_str(), "wb0");
  sam_hdr_write(fout, h_);
  
  int64_t record_length = 0;
  for (int i = 0; i < input.size; i++) {
    int tmp_l = bam_write1(fout->fp.bgzf, input.bams[i]);
    record_length += tmp_l;
    bam_destroy1(input.bams[i]);
  }
  free(input.bams);

  sam_close(fout);
  
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished BamWrite()";
}
