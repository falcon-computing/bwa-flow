#include <iostream>
#include <sstream>
#include <iomanip>

#include "BamReadStage.h"

BamRecord BamReadStage::compute(int const & id) {
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started BamRead()";
  std::stringstream ss;
  ss << temp_dir_ << "/part-" << std::setw(6) 
     << std::setfill('0') << id << ".bam";
  samFile * file_pointer = hts_open(ss.str().c_str(), "r");
  //samFile * file_pointer = hts_open("/genome/disk2/tianj/HG001.4.bam", "rb");
  //DLOG(INFO) << "check fp " << file_pointer->lineno;
  if (file_pointer == NULL) {
    throw("file not exist");
  }
  bam_hdr_t * bamHdr = sam_hdr_read(file_pointer);
  //DLOG(INFO) << "check header" << bamHdr->l_text;
  int align_size = 100000;
  bam1_t** aligns = (bam1_t**)malloc(align_size*sizeof(bam1_t*));
  int i = 0;  
  while (true) {
    bam1_t* align = bam_init1();
    int flag = sam_read1(file_pointer, bamHdr, align);
    if (flag < 0){
      //DLOG(INFO) << "file_path " << ss.str() ;
      //DLOG(INFO) << "flag value " << flag << " i: "<< i;
      bam_destroy1(align);
      sam_close(file_pointer);
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
  //DLOG(INFO) << "bucket "<< id << " i: " << i;
  output.id = id; 
  output.size = i;
  output.bams = aligns;
  
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished BamRead()";
  return output;
}
