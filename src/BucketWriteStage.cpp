#include "BucketWriteStage.h"
#include "util.h"

int BucketWriteStage::get_bucket_id(bam1_t* read) {
  int32_t contig_id = read->core.tid;
  int32_t read_pos = read->core.pos;
//DLOG(INFO) << "read_pos " << contig_id;
//DLOG(INFO) << "contig_id " << contig_id;
  int64_t acc_pos = _accumulate_length[contig_id] + read_pos;
//DLOG(INFO) << "acc_pos " << acc_pos;
  return (acc_pos-1)/_bucket_size;
}

int BucketWriteStage::compute(BamsRecord const & input) {
  uint64_t start = getUs();
  DLOG(INFO) << "Started BucketWrite()";
  std::unordered_map<int32_t, std::vector<bam1_t*> > buckets;
  for (int k = 0; k < input.records_list->size(); k++) {
    int batch_num = input.records_list[0][k].batch_num;
//DLOG(INFO) << "batch_num " << batch_num;
    for (int i = 0; i < batch_num; i++) {
      for (int j = 0; j < input.records_list[0][k].seqs[i].bams->l; j++) {
        bam1_t* tmp = input.records_list[0][k].seqs[i].bams->bams[j];
        int bucket_id = get_bucket_id(tmp);
        if (buckets.count(bucket_id) != 1) {
          std::vector<bam1_t*> tmp_vec;
          buckets[bucket_id] = tmp_vec;
        }
        buckets[bucket_id].push_back(tmp);
      }
    }
  }
  for (int i = 0; i < _buckets.size(); i++) {
    if (buckets.count(i)) {
      _buckets[i][0].writeFile(buckets[i]);
    }
    for (int j = 0; j < buckets[i].size(); j++) {
      bam_destroy1(buckets[i][j]);
    }
  }
  for (int k = 0; k < input.records_list->size(); k++) {
    int batch_num = input.records_list[0][k].batch_num;
    for (int i = 0; i < batch_num; i++) {
      free(input.records_list[0][k].seqs[i].bams->bams);
      free(input.records_list[0][k].seqs[i].bams);
      input.records_list[0][k].seqs[i].bams = NULL;
    }
  free(input.records_list[0][k].seqs);
  }
//  free(input.bam_buffer);
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