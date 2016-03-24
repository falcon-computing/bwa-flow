#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char** argv) {

  int n = 8;
  int length = 8;

  if (argc > 1) {
    n = atoi(argv[1]);
  }
  if (argc > 2) {
    length = atoi(argv[2]);
  }

  srand(time(NULL));

  ofstream fout_data("input.txt");
  ofstream fout_result("result.txt");

  for (int i=0; i<n; i++) {
    double norm = 0;
    for (int i=0; i<length; i++) {
      double val = (double)rand()/RAND_MAX;
      fout_data << val;
      if (i<length-1) fout_data << " ";
      norm += val*val;
    }
    fout_data << endl;
    fout_result << norm << endl;
  }
  
  return 0;
}
