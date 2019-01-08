#include <iostream>
#include <sstream>
#include <iomanip>

#include "BamReadStage.h"

BamRecord BamReadStage::compute(int const & id) {

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started BamRead()";

  std::stringstream ss;
  ss << bam_dir_ << "/part-" << std::setw(6) 
     << std::setfill('0') << id << ".bam";

  samFile * fp = hts_open(ss.str().c_str(), "r");

  if (!fp) {
    throw std::runtime_error("bucket bam not exist");
  }

  bam_hdr_t* h = sam_hdr_read(fp);
  if (!h) throw std::runtime_error("failed to read bam header");

  int align_size = 100000;
  bam1_t** aligns = (bam1_t**)malloc(align_size*sizeof(bam1_t*));

  int i = 0;  
  while (true) {
    bam1_t* align = bam_init1();
    if (sam_read1(fp, h, align) < 0) {
      bam_destroy1(align);
      break;
    }   

    if (i >= align_size) {
      align_size *= 2;
      aligns = (bam1_t**)realloc(aligns, align_size*sizeof(bam1_t*));
    }   
    aligns[i] = align;
    i++;
  }
  bam_hdr_destroy(h);
  sam_close(fp);

  BamRecord output;
  output.id = id; 
  output.size = i;
  output.bams = aligns;
  
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished BamRead()";
  return output;
}
