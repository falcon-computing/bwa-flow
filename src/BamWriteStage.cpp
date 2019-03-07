#include <boost/filesystem.hpp>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "BamWriteStage.h"
#include "config.h"
#include "util.h"

namespace bfs = boost::filesystem;

BamWriteStage::BamWriteStage(
    int num_parts,
    std::string bam_dir, 
    std::string output_path, 
    bam_hdr_t* h,
    int n): kestrelFlow::MapStage<
      BamRecord, int, COMPUTE_DEPTH, COMPUTE_DEPTH>(n), 
      num_parts_(num_parts),
      bam_dir_(bam_dir), 
      output_path_(output_path),
      h_(h)
{
  std::stringstream ss;
  ss << bam_dir_ << "/header";
  samFile* fout = hts_open(ss.str().c_str(), "wb");
  sam_hdr_write(fout, h_);
  sam_close(fout);

  // use boost to resize the file to remove EOF marker
  bfs::resize_file(ss.str(), 
      bfs::file_size(ss.str()) - 28);
} 

int BamWriteStage::compute(BamRecord const & input) {
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started BamWrite()";

  int id = input.id;
  std::stringstream ss;
  ss << bam_dir_ << "/part-" 
     << std::setw(6) << std::setfill('0') << id;
  
  samFile * fout = hts_open(ss.str().c_str(), "wb");

  // don't write headers for each part bam
  //sam_hdr_write(fout, h_);
  
  uint64_t start_ts = getUs();

  bgzf_write(fout->fp.bgzf, 
      input.fbuf->get_data(), 
      input.fbuf->get_size());

  delete input.fbuf;
  DLOG(INFO) << "Wrtting " << input.fbuf->get_size() / 1024 << " kB took " 
    << getUs() - start_ts << " us";

  sam_close(fout);
  
  // use boost to resize the file to remove EOF marker
  // leave the EOF marker for the last bucket part,
  // since it acts as the EOF for the entire file
  if (id != num_parts_ - 1) {
    bfs::resize_file(ss.str(), 
        bfs::file_size(ss.str()) - 28);
  }
  
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished BamWrite()";
}

BamWriteStage::~BamWriteStage() {
  // TODO: check if the command will be too long
  // if so, we need to cat separately to append
  char path_buf[2000];
  realpath(output_path_.c_str(), path_buf);
  
  std::string ab_output = std::string(path_buf);

  std::stringstream ss;
  ss << "cd " << bam_dir_ << " && ";
  ss << "cat " << "./header ";
  for (int i = 0; i < num_parts_; ++i) {
    ss << "./part-" 
       << std::setw(6) << std::setfill('0') << i
       << " ";
  }
  ss << "> " << ab_output;

  uint64_t start_ts = getUs();

  int ret = system(ss.str().c_str());

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "concat BAM files took " 
    << getUs() - start_ts << " us";
}
