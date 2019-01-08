#ifndef BUCKETSORT
#define BUCKETSORT

#include <cstring>
#include <string.h>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sys/stat.h>

#include <boost/atomic.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread.hpp>

#include "Pipeline.h"
#include "config.h"

//locakable class
class bucketFile :
  public boost::basic_lockable_adapter<boost::mutex>{
  private:
    samFile * fout_;
    ktp_aux_t * aux_;
    int32_t id_;
    char* file_path_;
    char* mode_;
    void writeFileHeader();
  public:
    int32_t get_id() {
      return id_;
    }
    bucketFile(ktp_aux_t* aux, int32_t id, const char* file_path,
              const char* mode): aux_(aux), id_(id) {
      file_path_ = strdup(file_path);
      mode_ = strdup(mode);
      fout_ = sam_open(file_path_, mode_);
      writeFileHeader();
    }
    ~bucketFile() {
      sam_close(fout_);
      free(file_path_);
      free(mode_);
    }
    void writeFile(std::vector<bam1_t*> vec);
};

class BucketSortStage :
  public kestrelFlow::MapStage<BamsRecord, int, COMPUTE_DEPTH, 0> {
  public:
    BucketSortStage(ktp_aux_t* aux, std::string out_dir, int num_buckets = 1, int n = 1):
      kestrelFlow::MapStage<BamsRecord, int, COMPUTE_DEPTH, 0>(n), aux_(aux) {
        accumulate_length_.push_back(0);
        int64_t acc_len = 0;
        for (int i = 0; i < aux_->h->n_targets; i++) {
          acc_len += aux_->h->target_len[i];
          accumulate_length_.push_back(acc_len);
        }
        const char *modes[] = {"wb", "wb0", "w"};
        for (int i = 0; i < num_buckets; i++) {
          //boost::any var = this->getConst("sam_dir");
          //std::string out_dir = boost::any_cast<std::string>(var);
          std::stringstream ss; 
          ss << out_dir << "/part-" << std::setw(6) << std::setfill('0') << i << ".bam";
          bucketFile* tmp_bucket = new bucketFile(aux_, i, ss.str().c_str(), modes[FLAGS_output_flag]);
          buckets_[i] = tmp_bucket;
        }
        std::stringstream ss;
        ss << out_dir << "/unmap.bam";
        star_read_ = new bucketFile(aux_, -1, ss.str().c_str(), modes[FLAGS_output_flag]);
        bucket_size_ = accumulate_length_[aux_->h->n_targets]/num_buckets;
        if (bucket_size_ == 0) {
          throw "bucket_size_ is 0";
        }
        if (accumulate_length_[aux_->h->n_targets]%num_buckets != 0) {
          bucket_size_ += 1;
        }
        std::stringstream interval_folder_path;
        interval_folder_path << out_dir << "/intervals";
        mkdir(interval_folder_path.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        std::ofstream interval_file;
        int bucket_per_thread = num_buckets/n;
        if (num_buckets%n != 0) {
          bucket_per_thread += 1;
        }
        int64_t thread_size = bucket_per_thread * bucket_size_;
        //int64_t thread_size = accumulate_length_[aux_->h->n_targets]/n;
        //if (accumulate_length_[aux_->h->n_targets]%n != 0) {
        //  thread_size +=1;
        //}
        int contig_start_pos = 0;
        int contig_id = 0;
        for (int i = 0; i < n; i++) {
          std::stringstream interval_file_path;
          interval_file_path << interval_folder_path.str().c_str();
          interval_file_path << "/interval_" << std::to_string(i) <<
            "_" << std::setw(6) << std::setfill('0') << i*bucket_per_thread << 
            "_" << std::setw(6) << std::setfill('0') << 
            std::min((i+1)*bucket_per_thread -1, num_buckets - 1) <<
            ".bed";
          interval_file.open(interval_file_path.str().c_str());
          int64_t end = contig_start_pos + thread_size;
          int margin_size = 1000;
          while (end > aux_->h->target_len[contig_id]) {
            // give some margin area to the interval
            int front_margin = 0;
            if (contig_start_pos >= margin_size) {
              front_margin = margin_size;
            }
            else {
              front_margin = contig_start_pos;
            }
            interval_file << aux_->h->target_name[contig_id] << "\t" << contig_start_pos - front_margin
              << "\t" << aux_->h->target_len[contig_id] << "\t" << i << "\n";
            end = end - aux_->h->target_len[contig_id];
            contig_start_pos = 0;
            contig_id += 1;
            if (contig_id >= aux_->h->n_targets) {
              DLOG(INFO) << "unexpected contig id exceeded.";
              break;
            }
          }
          if (contig_id >= aux_->h->n_targets) {
            break;
          }
          int end_margin = 0;
          if (end <= aux_->h->target_len[contig_id] - margin_size) {
            end_margin = margin_size;
          }
          else {
            end_margin = aux_->h->target_len[contig_id] - end;
          }
          interval_file << aux_->h->target_name[contig_id] << "\t" << contig_start_pos << "\t"
            << end + end_margin << "\t" << i << "\n";
          contig_start_pos = end;
          interval_file.close();
        }
        //interval_file.open(interval_file_path.str().c_str());
        //int contig_start_pos = 0;
        //int contig_id = 0;
        //for (int i = 0; i < num_buckets && contig_id < aux_->h->n_targets; i++) {
        //  int end = contig_start_pos + bucket_size_;
        //  while (end > aux_->h->target_len[contig_id]) {
        //    interval_file << aux_->h->target_name[contig_id] << "\t" << contig_start_pos 
        //      << "\t" << aux_->h->target_len[contig_id] << "\t" << i << "\n";
        //    end = end - aux_->h->target_len[contig_id];
        //    contig_start_pos = 0;
        //    contig_id += 1;
        //    if (contig_id >= aux_->h->n_targets) {
        //      DLOG(INFO) << "unexpected contig id exceeded.";
        //      break;
        //    }
        //  }
        //  if (contig_id >= aux_->h->n_targets) {
        //    break;
        //  }
        //  interval_file << aux_->h->target_name[contig_id] << "\t" << contig_start_pos << "\t"
        //    << end << "\t" << i << "\n";
        //  contig_start_pos = end;
        //}
        //interval_file.close();
      }
    ~BucketSortStage() {
//        for (auto it = buckets_.begin(); it != buckets_.end(); ++it) {
//          delete it->second;
//        }
    }
    int compute(BamsRecord const & input);
    void closeBuckets();
  private:
    ktp_aux_t* aux_;
    std::unordered_map<int32_t, bucketFile*> buckets_;
    int64_t bucket_size_;
    std::vector<int64_t> accumulate_length_;
    int get_bucket_id(bam1_t* read);
    bucketFile* star_read_;
};

#endif
