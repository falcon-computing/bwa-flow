#include "BucketWriteStage.h"

int BucketWriteStage::compute(SeqsRecord const & input) {
  std::unordered_map<int32_t, std::vector<bam1_t*> > buckets;
  int batch_num = input.batch_num;
  for (int i = 0; i < batch_num; i++) {
    bam1_t* tmp = input.seqs->bams->bams[i];
    if (buckets.count(tmp->core.tid) != 1) {
      std::vector<bam1_t*> tmp_vec;
      buckets[tmp->core.tid] = tmp_vec;
    }
    buckets[tmp->core.tid].push_back(tmp);
  }
  for (int i = 0; i < _aux->h->n_targets; i++) {
    _buckets[i][0].writeFile(&buckets[i][0], buckets[i].size());
  }
  return 0;
}

void bucketFile::writeFileHeader() {
  int status = sam_hdr_write(_fout, _aux->h);
  if (status) {
    ;
  }
  return;
}

void bucketFile::writeFile(bam1_t** reads, int batch_num) {
  boost::lock_guard<bucketFile> guard(*this);
  for (int i = 0; i < batch_num; i++) {
    sam_write1(_fout, _aux->h, reads[i]);
  }
  return;
}
