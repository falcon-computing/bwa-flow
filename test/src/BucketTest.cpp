#include "BucketSortStage.h"
#include "Pipeline.h"
#include "TestCommon.h"
#include "Pipeline.h"
#include "config.h"

TEST_F(BucketTest, intervalTest) {
  try {
    ktp_aux_t aux;
    bam_hdr_t headr;
    aux.h = &headr;
    aux.h->n_targets = 5;
    uint32_t len_tmp[5] = {1,2,3,2,1};
    aux.h->target_len = len_tmp;
    char * cc[] = {"0", "1", "2", "3", "4"};
    aux.h->target_name = cc;

    BucketSortStage stage(&aux, 3);
    ASSERT_EQ(stage.get_bucket_size(), 3);

    std::vector<std::vector<int64_t>> vec_vec = stage.get_intervals(0, 1);
    ASSERT_EQ(vec_vec.size(), 1);
    ASSERT_EQ(vec_vec[0].size(), 3);
    ASSERT_EQ(vec_vec[0][0], 0);
    ASSERT_EQ(vec_vec[0][1], 0);
    ASSERT_EQ(vec_vec[0][2], 1);
    vec_vec = stage.get_intervals(3, 5);
    ASSERT_EQ(vec_vec.size(), 1);
    ASSERT_EQ(vec_vec[0][0], 2);
    ASSERT_EQ(vec_vec[0][1], 0);
    ASSERT_EQ(vec_vec[0][2], 2);
    vec_vec = stage.get_intervals(1, 9);
    ASSERT_EQ(vec_vec.size(), 4);
    ASSERT_EQ(vec_vec[0][0], 1);
    ASSERT_EQ(vec_vec[0][1], 0);
    ASSERT_EQ(vec_vec[0][2], 2);
    ASSERT_EQ(vec_vec[1][0], 2);
    ASSERT_EQ(vec_vec[1][1], 0);
    ASSERT_EQ(vec_vec[1][2], 3);
    ASSERT_EQ(vec_vec[2][0], 3);
    ASSERT_EQ(vec_vec[2][1], 0);
    ASSERT_EQ(vec_vec[2][2], 2);
    ASSERT_EQ(vec_vec[3][0], 4);
    ASSERT_EQ(vec_vec[3][1], 0);
    ASSERT_EQ(vec_vec[3][2], 1);
    int index = stage.bucket_id_calculate(0, 0);
    ASSERT_EQ(index, 0);
    index = stage.bucket_id_calculate(4, 0);
    ASSERT_EQ(index, 2);
    index = stage.bucket_id_calculate(2, 2);
    ASSERT_EQ(index, 1);
  }
  catch (...) {
    FAIL() << "should not catch exception";
  }

  try {
    ktp_aux_t aux;
    bam_hdr_t headr;
    aux.h = &headr;
    aux.h->n_targets = 5;
    uint32_t len_tmp[5] = {4,3,2,1,1};
    aux.h->target_len = len_tmp;
    char * cc[] = {"0", "1", "2", "3", "4"};
    aux.h->target_name = cc;

    BucketSortStage stage(&aux, 10);
    ASSERT_EQ(stage.get_bucket_size(), 2);

    std::vector<std::vector<int64_t>> vec_vec = stage.get_intervals(0, 1);
    ASSERT_EQ(vec_vec.size(), 1);
    ASSERT_EQ(vec_vec[0].size(), 3);
    ASSERT_EQ(vec_vec[0][0], 0);
    ASSERT_EQ(vec_vec[0][1], 0);
    ASSERT_EQ(vec_vec[0][2], 1);
    vec_vec = stage.get_intervals(3, 5);
    ASSERT_EQ(vec_vec.size(), 2);
    ASSERT_EQ(vec_vec[0][0], 0);
    ASSERT_EQ(vec_vec[0][1], 3);
    ASSERT_EQ(vec_vec[0][2], 4);
    ASSERT_EQ(vec_vec[1][0], 1);
    ASSERT_EQ(vec_vec[1][1], 0);
    ASSERT_EQ(vec_vec[1][2], 1);
    vec_vec = stage.get_intervals(1, 8);
    ASSERT_EQ(vec_vec.size(), 3);
    ASSERT_EQ(vec_vec[0][0], 0);
    ASSERT_EQ(vec_vec[0][1], 1);
    ASSERT_EQ(vec_vec[0][2], 4);
    ASSERT_EQ(vec_vec[1][0], 1);
    ASSERT_EQ(vec_vec[1][1], 0);
    ASSERT_EQ(vec_vec[1][2], 3);
    ASSERT_EQ(vec_vec[2][0], 2);
    ASSERT_EQ(vec_vec[2][1], 0);
    ASSERT_EQ(vec_vec[2][2], 1);
    int index = stage.bucket_id_calculate(0, 0);
    ASSERT_EQ(index, 0);
    index = stage.bucket_id_calculate(0, 1);
    ASSERT_EQ(index, 0);
    index = stage.bucket_id_calculate(3, 0);
    ASSERT_EQ(index, 8);
    index = stage.bucket_id_calculate(1, 2);
    ASSERT_EQ(index, 5);
  }
  catch (...) {
    FAIL() << "should not catch exception";
  }

  try {
    ktp_aux_t aux;
    bam_hdr_t headr;
    aux.h = &headr;
    aux.h->n_targets = 5;
    uint32_t len_tmp[5] = {4,3,2,1,1};
    aux.h->target_len = len_tmp;
    char * cc[] = {"0", "1", "2", "3", "4"};
    aux.h->target_name = cc;

    BucketSortStage stage(&aux, 20);
    ASSERT_EQ(stage.get_bucket_size(), 1);

    int index = stage.bucket_id_calculate(0, 0);
    ASSERT_EQ(index, 0);
    index = stage.bucket_id_calculate(0, 1);
    ASSERT_EQ(index, 1);
    index = stage.bucket_id_calculate(3, 0);
    ASSERT_EQ(index, 9);
    index = stage.bucket_id_calculate(1, 2);
    ASSERT_EQ(index, 6);
  }
  catch (...) {
    FAIL() << "should not catch exception";
  }

}
