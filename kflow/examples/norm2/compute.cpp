#include <iostream>
#include <fstream>
#include <vector>
#include "kflow/Pipeline.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"
#include "kflow/MapStage.h"

using namespace kestrelFlow;

template <typename T>
class Record {
  public:
    Record(int _id): id_(_id) {}
    Record(int _id, T _data): id_(_id), data(_data) {}
    int id() { return id_; }
    T data;
  private:
    int id_;
};

class Load : 
  public SourceStage<Record<std::vector<double> >*>
{
public:
  Load(): SourceStage<Record<std::vector<double> >*>() {;}

  void compute() {

    try {
      boost::any var = this->getConst("input_fname");
      std::string fname = boost::any_cast<std::string>(var);

      std::ifstream fin(fname);

      std::string line;
      int idx = 0;
      while (std::getline(fin, line)) {
        std::stringstream ss(line);

        Record<std::vector<double> >* record = 
          new Record<std::vector<double> >(idx);

        double val;
        while (ss >> val) {
          record->data.push_back(val);
        }
        pushOutput(record);
        idx ++;
      }
    }
    catch (paramError &e) {
      LOG(ERROR) << e.what();
      return;
    }
    catch (boost::bad_any_cast &) {
      LOG(ERROR) << "type of input_fname mismatch";
      return;
    }
  }
};

template <int BSIZE = 64>
class Store : 
  public SinkStage<Record<double>*, BSIZE>
{
  public:
    Store(): SinkStage<Record<double>*, BSIZE>() {}

    void compute() {
      try {
        boost::any var = this->getConst("output_fname");
        std::string fname = boost::any_cast<std::string>(var);

        std::ofstream fout(fname);

        int counter = 0;
        std::vector<double> reorder_buffer(BSIZE);
        while (true) {

          Record<double>* input;
          bool ready = this->getInput(input);

          while (!this->isFinal() && !ready) {
            boost::this_thread::sleep_for(boost::chrono::microseconds(100));
            ready = this->getInput(input);
          }
          if (!ready) { 
            // this means isFinal() is true and input queue is empty
            break; 
          }
          reorder_buffer[input->id() % BSIZE] = input->data;
          counter ++;
          if (counter == BSIZE) {
            // write file
            for (int i=0; i<BSIZE; i++) {
              fout << reorder_buffer[i] << "\n";
            }
            // reset counter
            counter = 0;
          }
        }
        // write the last batch
        if (counter > 0) {
          for (int i=0; i<counter; i++) {
            fout << reorder_buffer[i] << "\n";
          }
        }
      }
      catch (paramError &e) {
        LOG(ERROR) << e.what();
        return;
      }
      catch (boost::bad_any_cast &) {
        LOG(ERROR) << "type of input_fname mismatch";
        return;
      }
    }
};

class Norm :
  public MapStage<
    Record<std::vector<double> >*, 
    Record<double>*>
{
public:
  Norm(bool is_dyn=true, int n_workers=1): 
    MapStage<
      Record<std::vector<double> >*, 
      Record<double>* 
    >(n_workers, is_dyn) {}

  Record<double>* compute(
      Record<std::vector<double> >* const & input) 
  {
    double norm = 0;
    for (int i=0; i<input->data.size(); i++) {
      double val = input->data[i];
      norm += val*val;
    }
    return new Record<double>(input->id(), norm);
  }
};


int main(int argc, char** argv) {

  FLAGS_logtostderr = 1;
  google::InitGoogleLogging(argv[0]);

  Pipeline norm(3, 0);

  std::string input_fname = "input.txt";
  std::string output_fname = "output.txt";
  norm.addConst("input_fname", input_fname);
  norm.addConst("output_fname", output_fname);

  Load  stage0;
  Norm  stage1(false, 2);
  Store<> stage2;

  norm.addStage(0, &stage0);
  norm.addStage(1, &stage1);
  norm.addStage(2, &stage2);
  norm.start();

  // gracefully end the pipeline
  norm.wait();

  return 0;
}
