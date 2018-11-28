#include "BucketWriteStage.h"
#include "util.h"

int BucketWriteStage::compute(BamsRecord const & input) {
  uint64_t start = getUs();
  DLOG(INFO) << "Started BucketWrite()";
  std::unordered_map<int32_t, std::vector<bam1_t*> > buckets;
  for (int k = 0; k < input.records_list->size(); k++) {
    int batch_num = input.records_list[0][k].batch_num;
    for (int i = 0; i < batch_num; i++) {
      for (int j = 0; j < input.records_list[0][k].seqs[i].bams->l; j++) {
        bam1_t* tmp = input.records_list[0][k].seqs[i].bams->bams[j];
        if (buckets.count(tmp->core.tid) != 1) {
          std::vector<bam1_t*> tmp_vec;
          buckets[tmp->core.tid] = tmp_vec;
        }
        buckets[tmp->core.tid].push_back(tmp);
      }
    }
  }
  for (int i = 0; i < _aux->h->n_targets; i++) {
    _buckets[i][0].writeFile(buckets[i]);
  }
  uint64_t end = getUs();
  DLOG(INFO) << "Finished BucketWrite() in " << start-end << " us.";
  return 0;
}

void bucketFile::writeFileHeader() {
  int status = sam_hdr_write(_fout, _aux->h);
  if (status) {
    ;
  }
  return;
}

void bucketFile::writeFile(std::vector<bam1_t*> vec) {
  boost::lock_guard<bucketFile> guard(*this);
  _fout = sam_open(_file_path, _mode);
  for (int i = 0; i < vec.size(); i++) {
    sam_write1(_fout, _aux->h, vec[i]);
  }
  return;
}
