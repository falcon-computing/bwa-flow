#ifndef SWCLIENT_H
#define SWCLIENT_H

#include "blaze/Client.h"
#include <glog/logging.h>

class SWClient : public blaze::Client {
 public:
  SWClient(): blaze::Client("SmithWaterman", 2, 1) {;}

  void compute() {
    throw std::runtime_error("No CPU implementation");
  }
};

#endif
