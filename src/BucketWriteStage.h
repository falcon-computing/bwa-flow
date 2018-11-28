#ifndef BUCKETWRITE
#define BUCKETWRITE

#include "Pipeline.h"
#include <cstring>
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
    void writeFileHeader();
  public:
    int32_t get_id() {
      return _id;
    }
    bucketFile(ktp_aux_t* aux, int32_t id, const char* file_path,
              const char* mode): _aux(aux), _id(id) {
      _fout = sam_open(file_path, mode);
      writeFileHeader();
    } 
    ~bucketFile() {
      sam_close(_fout);
    }
    void writeFile(std::vector<bam1_t*> vec);
};

class BucketWriteStage :
  public kestrelFlow::MapStage<SeqsRecord, int, COMPUTE_DEPTH, 0> {
  public:
    BucketWriteStage(ktp_aux_t* aux, std::string out_dir, int n = 1):
      kestrelFlow::MapStage<SeqsRecord, int, COMPUTE_DEPTH, 0>(n, false), _aux(aux) {
        // initialize buckets
//DLOG(INFO) << "constructor BW";
        for (int i = 0; i < _aux->h->n_targets; i++) {
//DLOG(INFO) << "constructing " << i << "  buckets";
          //boost::any var = this->getConst("sam_dir");
          //std::string out_dir = boost::any_cast<std::string>(var);
          std::stringstream ss;    
          ss << out_dir << "/contig-" << std::setw(3) << std::setfill('0') << i << ".bam";
          const char *modes[] = {"wb", "wb0", "w"};
          bucketFile* tmp_bucket = new bucketFile(_aux, i, ss.str().c_str(), modes[FLAGS_output_flag]);
          _buckets[i] = tmp_bucket;
        }
//DLOG(INFO) << "id of bucket 2 " << _buckets[2]->get_id();
      }
    ~BucketWriteStage() {
        for (auto it = _buckets.begin(); it != _buckets.end(); ++it) {
          delete it->second;
        }
    }
    int compute(SeqsRecord const & input);
  private:
    ktp_aux_t* _aux;
    std::unordered_map<int32_t, bucketFile*> _buckets;
};

#endif
