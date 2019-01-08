#include "BucketSortStage.h"
#include "util.h"

#define UNMAP_FLAG 4
#define DUP_FLAG 1024

int BucketSortStage::get_bucket_id(bam1_t* read) {
  if (read->core.tid == -1) {

    return num_buckets_;
  }
  int32_t contig_id = read->core.tid;
  int32_t read_pos = read->core.pos;
  int64_t acc_pos = accumulate_length_[contig_id] + read_pos;
  return (acc_pos-1)/bucket_size_;
}

void BucketSortStage::closeBuckets() {
  for (auto it = buckets_.begin(); it != buckets_.end(); ++it) {
    delete it->second;
  }
  //delete star_read_;
}

int BucketSortStage::compute(BamsRecord const & input) {
  uint64_t start = getUs();
  DLOG(INFO) << "Started BucketWrite()";
  std::unordered_map<int32_t, std::vector<bam1_t*> > buckets;

  for (int k = 0; k < input.records_list->size(); k++) {
    int batch_num = input.records_list[0][k].batch_num;
    for (int i = 0; i < batch_num; i++) {
      for (int j = 0; j < input.records_list[0][k].seqs[i].bams->l; j++) {
        bam1_t* tmp = input.records_list[0][k].seqs[i].bams->bams[j];
        int bucket_id = get_bucket_id(tmp);

        buckets[bucket_id].push_back(tmp);
      }
    }
  }
  for (int i = 0; i < buckets_.size(); i++) {
    if (buckets.count(i)) {
      buckets_[i][0].writeFile(buckets[i]);
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
  DLOG(INFO) << "Finished BucketWrite() in " << getUs() - start << " us.";
  return 0;
}

void bucketFile::writeFileHeader() {
  int status = sam_hdr_write(fout_, aux_->h);
  if (status) {
    ;
  }
  return;
}

void bucketFile::writeFile(std::vector<bam1_t*> vec) {
  boost::lock_guard<bucketFile> guard(*this);
  for (int i = 0; i < vec.size(); i++) {
    if (!(FLAGS_remove_duplicates && (vec[i]->core.flag & DUP_FLAG))) {
      sam_write1(fout_, aux_->h, vec[i]);
    }
  }
}
