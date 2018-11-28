#include "BucketWriteStage.h"
#include "util.h"

int BucketWriteStage::compute(SeqsRecord const & input) {
  uint64_t start = getUs();
  DLOG(INFO) << "Started BucketWrite()";
  std::unordered_map<int32_t, std::vector<bam1_t*> > buckets;
  int batch_num = input.batch_num;
  for (int i = 0; i < batch_num; i++) {
    for (int j = 0; j < input.seqs[i].bams->l; j++) {
      bam1_t* tmp = input.seqs[i].bams->bams[j];
      if (buckets.count(tmp->core.tid) != 1) {
        std::vector<bam1_t*> tmp_vec;
        buckets[tmp->core.tid] = tmp_vec;
      }
      buckets[tmp->core.tid].push_back(tmp);
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
  for (int i = 0; i < vec.size(); i++) {
    sam_write1(_fout, _aux->h, vec[i]);
  }
  return;
}
