#ifndef COMMON_H
#define COMMON_H

#include <boost/asio.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <cstdint>
#include <glog/logging.h>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace kestrelFlow 
{

class Pipeline;
class StageBase;
class QueueBase;

// base class of a data record
class RecordBase {
  public: 
    RecordBase(int _s = 0, int _n = 0): 
      size_(_s), n_elements_(_n) {}

    virtual ~RecordBase() {}

  protected:
    int size_;
    int n_elements_;
};

// custom exceptions
class paramError : public std::runtime_error {
public:
  explicit paramError(const std::string& what_arg):
    std::runtime_error(what_arg) {;}
};

class fileError : public std::runtime_error {
public:
  explicit fileError(const std::string& what_arg):
    std::runtime_error(what_arg) {;}
};

class internalError : public std::runtime_error {
public:
  explicit internalError(const std::string& what_arg):
    std::runtime_error(what_arg) {;}
};

} // namespace kestrelFlow

#endif
