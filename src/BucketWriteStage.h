#ifndef BUCKETWRITE
#define BUCKETWRITE

#include "Pipeline.h"
#include <cstring>
#include <string.h>
#include <unordered_map>
#include "config.h"

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

class BucketWriteStage :
  public kestrelFlow::MapStage<BamsRecord, int, COMPUTE_DEPTH, 0> {
  public:
    BucketWriteStage(ktp_aux_t* aux, std::string out_dir, int num_buckets = 1, int n = 1):
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
//DLOG(INFO) << "id of bucket 2 " << _buckets[2]->get_id();
      }
    ~BucketWriteStage() {
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
