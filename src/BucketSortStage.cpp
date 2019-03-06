#include "BucketSortStage.h"
#include "util.h"

#include <string>
#include <vector>
#include <sstream>

#define UNMAP_FLAG 4
#define DUP_FLAG 1024

std::vector<std::vector<int64_t>> BucketSortStage::get_intervals(int64_t start, int64_t end) {
  if (start < 0 || end < 0) throw("pos less than 0");
  if (start > end) throw("start pos larger than end");
  if (start > accumulate_length_[aux_->h->n_targets]) throw("start pos out of reference boundary");
  std::vector<std::vector<int64_t>> res;
  std::vector<int64_t> tmp;
  // calculate the start contig id;
  int contig_id = 0;
  while (contig_id < aux_->h->n_targets && 
        start >= accumulate_length_[contig_id + 1]) {
    contig_id ++;
  }
  while (contig_id < aux_->h->n_targets &&
        end > accumulate_length_[contig_id + 1]) {
    tmp.push_back(contig_id);
    tmp.push_back(start - accumulate_length_[contig_id]);
    tmp.push_back(accumulate_length_[contig_id + 1] - accumulate_length_[contig_id]);
    res.push_back(tmp);
    tmp.clear();
    start = accumulate_length_[contig_id + 1];
    contig_id ++;
  }
  if (contig_id < aux_->h->n_targets) {
    tmp.push_back(contig_id);
    tmp.push_back(start - accumulate_length_[contig_id]);
    tmp.push_back(end - accumulate_length_[contig_id]);
    res.push_back(tmp);
    tmp.clear();
  }
  return res;
}

int BucketSortStage::bucket_id_calculate(int32_t contig_id, int32_t read_pos) {
  int64_t acc_pos = accumulate_length_[contig_id] + read_pos;
  int large_bucket = (accumulate_length_[aux_->h->n_targets]%num_buckets_)?
                      (accumulate_length_[aux_->h->n_targets]%num_buckets_):
                      num_buckets_;
  int64_t limit = large_bucket * bucket_size_;
  return (acc_pos > limit)?
        (
          (bucket_size_ - 1)?
          (large_bucket + (acc_pos - limit)/(bucket_size_ - 1)):
          (large_bucket)
        ):
        (acc_pos/bucket_size_);
}

int BucketSortStage::get_bucket_id(bam1_t* read) {
  if (read->core.tid == -1) {

    return num_buckets_;
  }
  int32_t contig_id = read->core.tid;
  int32_t read_pos = read->core.pos;
  return bucket_id_calculate(contig_id, read_pos);
}

void BucketSortStage::closeBuckets() {
  for (auto it = buckets_.begin(); it != buckets_.end(); ++it) {
    delete it->second;
  }
  //delete star_read_;
}

BucketSortStage::BucketSortStage(ktp_aux_t* aux, int num_buckets):
  aux_(aux), num_buckets_(num_buckets) {
  accumulate_length_.push_back(0);
  int64_t acc_len = 0;
  for (int i = 0; i < aux_->h->n_targets; i++) {
    acc_len += aux_->h->target_len[i];
    accumulate_length_.push_back(acc_len);
  }
  bucket_size_ = (accumulate_length_[aux_->h->n_targets] + num_buckets - 1)/num_buckets;
}

BucketSortStage::BucketSortStage(
        ktp_aux_t* aux,
        std::string out_dir,
        int num_buckets,
        int n):
      kestrelFlow::MapStage<BamsRecord, int, COMPUTE_DEPTH, 0>(n),
      aux_(aux),
      num_buckets_(num_buckets) {
  accumulate_length_.push_back(0);
  int64_t acc_len = 0;
  for (int i = 0; i < aux_->h->n_targets; i++) {
    acc_len += aux_->h->target_len[i];
    accumulate_length_.push_back(acc_len);
  }
  const char *modes[] = {"wb", "wb0", "w"};
  int64_t contig_start_pos = 0;
  int contig_id = 0;
  bucket_size_ = (accumulate_length_[aux_->h->n_targets] + num_buckets - 1)/num_buckets;
  int large_bucket = accumulate_length_[aux_->h->n_targets]%num_buckets;
  for (int i = 0; i <= num_buckets; i++) {
    std::stringstream ss;
    ss << out_dir << "/part-" << std::setw(6) << std::setfill('0') << i << ".bam";
    buckets_[i] = new bucketFile(aux_, i, ss.str().c_str(), modes[FLAGS_output_flag]);

    if (i == num_buckets) break;

    std::stringstream ss2;
    ss2 << out_dir << "/part-" << std::setw(6) << std::setfill('0') << i << ".bed";
    std::ofstream intv_file(ss2.str().c_str());
    int64_t end = contig_start_pos + bucket_size_ - (int)(i >= large_bucket);
    std::vector<std::vector<int64_t>> intv_vec_vec = get_intervals(contig_start_pos, end);
    for (auto & intv_vec : intv_vec_vec) {
      intv_file << aux_->h->target_name[intv_vec[0]] << "\t"
                << intv_vec[1] << "\t"
                << intv_vec[2] << "\n";
    }
    contig_start_pos = end;
    intv_file.close();
  }
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
