#ifndef BUCKETSORT
#define BUCKETSORT

#include "Pipeline.h"
#include <cstring>
#include <string.h>
#include <unordered_map>
#include "config.h"

#include <iostream>
#include <fstream>
#include <sys/stat.h>

#include <boost/atomic.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread.hpp>

//locakable class
class bucketFile :
  public boost::basic_lockable_adapter<boost::mutex>{
  private:
    samFile * _fout;
    ktp_aux_t * _aux;
    int32_t _id;
    char* _file_path;
    char* _mode;
    void writeFileHeader();
  public:
    int32_t get_id() {
      return _id;
    }
    bucketFile(ktp_aux_t* aux, int32_t id, const char* file_path,
              const char* mode): _aux(aux), _id(id) {
      _file_path = strdup(file_path);
      _mode = strdup(mode);
      _fout = sam_open(_file_path, _mode);
      writeFileHeader();
    }
    ~bucketFile() {
      sam_close(_fout);
      free(_file_path);
      free(_mode);
    }
    void writeFile(std::vector<bam1_t*> vec);
};

class BucketSortStage :
  public kestrelFlow::MapStage<BamsRecord, int, COMPUTE_DEPTH, 0> {
  public:
    BucketSortStage(ktp_aux_t* aux, std::string out_dir, int num_buckets = 1, int n = 1):
      kestrelFlow::MapStage<BamsRecord, int, COMPUTE_DEPTH, 0>(n), _aux(aux) {
        _accumulate_length.push_back(0);
        int64_t acc_len = 0;
        for (int i = 0; i < _aux->h->n_targets; i++) {
          acc_len += _aux->h->target_len[i];
          _accumulate_length.push_back(acc_len);
        }
        for (int i = 0; i < num_buckets; i++) {
          //boost::any var = this->getConst("sam_dir");
          //std::string out_dir = boost::any_cast<std::string>(var);
          std::stringstream ss; 
          ss << out_dir << "/part-" << std::setw(6) << std::setfill('0') << i << ".bam";
          const char *modes[] = {"wb", "wb0", "w"};
          bucketFile* tmp_bucket = new bucketFile(_aux, i, ss.str().c_str(), modes[FLAGS_output_flag]);
          _buckets[i] = tmp_bucket;
        }
        _bucket_size = _accumulate_length[_aux->h->n_targets]/num_buckets;
        if (_bucket_size == 0) {
          throw "_bucket_size is 0";
        }
        if (_accumulate_length[_aux->h->n_targets]%num_buckets != 0) {
          _bucket_size += 1;
        }
        std::stringstream interval_folder_path;
        interval_folder_path << out_dir << "/intervals";
        mkdir(interval_folder_path.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        std::ofstream interval_file;
        int bucket_per_thread = num_buckets/n;
        if (num_buckets%n != 0) {
          bucket_per_thread += 1;
        }
        int64_t thread_size = bucket_per_thread * _bucket_size;
        //int64_t thread_size = _accumulate_length[_aux->h->n_targets]/n;
        //if (_accumulate_length[_aux->h->n_targets]%n != 0) {
        //  thread_size +=1;
        //}
        int contig_start_pos = 0;
        int contig_id = 0;
        for (int i = 0; i < n; i++) {
          std::stringstream interval_file_path;
          interval_file_path << interval_folder_path.str().c_str();
          interval_file_path << "/interval_" << std::to_string(i) << ".bed";
          interval_file.open(interval_file_path.str().c_str());
          int64_t end = contig_start_pos + thread_size;
          int margin_size = 1000;
          while (end > _aux->h->target_len[contig_id]) {
            // give some margin area to the interval
            int front_margin = 0;
            if (contig_start_pos >= margin_size) {
              front_margin = margin_size;
            }
            else {
              front_margin = contig_start_pos;
            }
            interval_file << _aux->h->target_name[contig_id] << "\t" << contig_start_pos - front_margin
              << "\t" << _aux->h->target_len[contig_id] << "\t" << i << "\n";
            end = end - _aux->h->target_len[contig_id];
            contig_start_pos = 0;
            contig_id += 1;
            if (contig_id >= _aux->h->n_targets) {
              DLOG(INFO) << "unexpected contig id exceeded.";
              break;
            }
          }
          if (contig_id >= _aux->h->n_targets) {
            break;
          }
          int end_margin = 0;
          if (end <= _aux->h->target_len[contig_id] - margin_size) {
            end_margin = margin_size;
          }
          else {
            end_margin = _aux->h->target_len[contig_id] - end;
          }
          interval_file << _aux->h->target_name[contig_id] << "\t" << contig_start_pos << "\t"
            << end + end_margin << "\t" << i << "\n";
          contig_start_pos = end;
          interval_file.close();
        }
        //interval_file.open(interval_file_path.str().c_str());
        //int contig_start_pos = 0;
        //int contig_id = 0;
        //for (int i = 0; i < num_buckets && contig_id < _aux->h->n_targets; i++) {
        //  int end = contig_start_pos + _bucket_size;
        //  while (end > _aux->h->target_len[contig_id]) {
        //    interval_file << _aux->h->target_name[contig_id] << "\t" << contig_start_pos 
        //      << "\t" << _aux->h->target_len[contig_id] << "\t" << i << "\n";
        //    end = end - _aux->h->target_len[contig_id];
        //    contig_start_pos = 0;
        //    contig_id += 1;
        //    if (contig_id >= _aux->h->n_targets) {
        //      DLOG(INFO) << "unexpected contig id exceeded.";
        //      break;
        //    }
        //  }
        //  if (contig_id >= _aux->h->n_targets) {
        //    break;
        //  }
        //  interval_file << _aux->h->target_name[contig_id] << "\t" << contig_start_pos << "\t"
        //    << end << "\t" << i << "\n";
        //  contig_start_pos = end;
        //}
        //interval_file.close();
      }
    ~BucketSortStage() {
        for (auto it = _buckets.begin(); it != _buckets.end(); ++it) {
          delete it->second;
        }
    }
    int compute(BamsRecord const & input);
  private:
    ktp_aux_t* _aux;
    std::unordered_map<int32_t, bucketFile*> _buckets;
    int64_t _bucket_size;
    std::vector<int64_t> _accumulate_length;
    int get_bucket_id(bam1_t* read);
};

#endif
