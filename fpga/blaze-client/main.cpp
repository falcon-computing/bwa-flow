#include <glog/logging.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "blaze/AccAgent.h"

#define FPGA_RET_PARAM_NUM 5

int main(int argc, char** argv) {
 
  char* dump_fname = argv[1];
  if (argc < 3) {
    printf("USAGE: %s <dump_file> <conf_file>\n", argv[0]);
    return 1;
  }

  try {
    std::string acc_id = "SmithWaterman";
    blaze::AccAgent agent(argv[2]);

    FILE* fin = fopen(dump_fname, "rb");
    if (!fin) {
      printf("Cannot find %s\n", dump_fname);
      return 1;
    }

    int data_size = 0;
    int task_num = 0;
    fread(&task_num, 1, sizeof(int), fin);
    fread(&data_size, 1, sizeof(int), fin);

    printf("reading %s, data_size=%d, task_num=%d\n", 
        dump_fname, data_size, task_num);

    int output_size  = FPGA_RET_PARAM_NUM*task_num;

    int* data_ptr   = new int[data_size];
    int* output_ptr = new int[output_size];

    // read from file
    fread(data_ptr, sizeof(int), data_size, fin);

    // start computation
    uint64_t start_ts = blaze::getUs();

    blaze::Task_ptr task = agent.createTask(acc_id);
    agent.writeInput(task, acc_id, data_ptr, 1, data_size, sizeof(int));
    agent.writeInput(task, acc_id, &task_num, 1, 1, sizeof(int));

    agent.readOutput(task, output_ptr, output_size*sizeof(int));

    printf("elapsed time for client: %dus\n", blaze::getUs()-start_ts);

    int* results = (int*) malloc(output_size*sizeof(int));
    fread(results, output_size, sizeof(int), fin);
    fclose(fin);

    // check results
    int correct = true;
    int err_count = 0;
    for (int k=0; k<output_size; k++) {
      if (output_ptr[k] != results[k]) {
        correct = false;
        //printf("%d != %d\n", output_ptr[k], results[k]);
        err_count++;
      }
    }
    if (correct) {
      printf("Results correct\n");
    }
    else {
      printf("%d/%d Results incorrect\n", err_count, output_size);
    }

    delete [] data_ptr;
    delete [] output_ptr;
  }
  catch (std::exception &e) {
    printf("%s\n", e.what());
    return -1;
  }

  return 0;
}
