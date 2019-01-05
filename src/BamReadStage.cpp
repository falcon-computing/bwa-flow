#include <iostream>
#include <sstream>
#include <iomanip>

#include "BamReadStage.h"

BamRecord BamReadStage::compute(int const & id) {
  std::stringstream ss;
  ss << temp_dir_ << "/part-" << std::setw(6) 
     << std::setfill('0') << id << ".bam";
  samFile * fp = hts_open(ss.str().c_str(), "r");
  if (fp == NULL) {
    throw("file not exist");
  }
  bam_hdr_t * bamHdr = sam_hdr_read(fp);
  int align_size = 100000;
  bam1_t** aligns = (bam1_t**)malloc(align_size*sizeof(bam1_t*));
  int i = 0;  
  while (true) {
    bam1_t* align = bam_init1();
    int flag = sam_read1(fp, bamHdr, align);
    if (flag == 0){ 
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
  BamRecord output;

  output.id = id; 
  output.size = i;
  output.bams = aligns;
  
  return output;
}
