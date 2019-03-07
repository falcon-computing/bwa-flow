#ifndef BUCKETSORT
#define BUCKETSORT

#include <cstring>
#include <string.h>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <gtest/gtest.h>
#include <gtest/gtest_prod.h>

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
    BucketSortStage( // for test only
        ktp_aux_t* aux, 
        int num_buckets);
  
    BucketSortStage(
        ktp_aux_t* aux, 
        std::string out_dir, 
        int num_buckets = 1, 
        int n = 1);
    
    ~BucketSortStage() {
//        for (auto it = buckets_.begin(); it != buckets_.end(); ++it) {
//          delete it->second;
//        }
    }
    int compute(BamsRecord const & input);
    void closeBuckets();
    int get_bucket_size() {
      return bucket_size_;
    }

  private:
    int num_buckets_;
    int64_t bucket_size_;
    ktp_aux_t* aux_;
    std::unordered_map<int32_t, bucketFile*> buckets_;
    std::vector<int64_t> accumulate_length_;
    FRIEND_TEST(BucketTest, intervalTest);
    
    int get_bucket_id(bam1_t* read);
    int bucket_id_calculate(int32_t contig_id, int32_t read_pos);
    std::vector<std::vector<int64_t>> get_intervals(int64_t start, int64_t end);
};

#endif
