#ifndef KFLOW_QUEUE_H
#define KFLOW_QUEUE_H

#include <boost/atomic.hpp>
#include <boost/lockfree/queue.hpp>
#include "Common.h"

namespace kestrelFlow 
{
class QueueBase {
 public:
   virtual bool empty() = 0;
};

template <typename U, int DEPTH = 64>
class Queue : public QueueBase {
 public:
  Queue<U, DEPTH>(int depth = DEPTH): num_elements_(0),
                                      capacity_(depth),
                                      data_queue_(depth) {}

  bool empty() {
    return (num_elements_.load() == 0);
  }

  bool almost_full() {
    return (num_elements_.load() >= capacity_ / 2);
  }

  int get_size() {
    return num_elements_.load();
  }

  int get_capacity() {
    return capacity_;
  }

  void pop(U &item) {
    while (!data_queue_.pop(item)) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(100));
    }
    num_elements_.fetch_sub(1);
  }

  void push(U item) {
    while (!data_queue_.push(item)) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(100));
    }
    num_elements_.fetch_add(1);
  }

  bool async_pop(U &item) {
    bool done = data_queue_.pop(item);
    if (done) num_elements_.fetch_sub(1);
    return done;
  }

  bool async_push(U item) {
    bool done = data_queue_.push(item);
    if (done) num_elements_.fetch_add(1);
    return done;
  }

 private:
  const int capacity_;
  mutable boost::atomic<int> num_elements_;
  //boost::lockfree::queue< U, boost::lockfree::capacity<DEPTH> > data_queue_;
  boost::lockfree::queue<U, boost::lockfree::fixed_sized<true> > data_queue_;
};

template <>
class Queue<void, 0> : public QueueBase {
 public:
   bool empty() { return true;}
};

typedef boost::shared_ptr<QueueBase> Queue_ptr;
const Queue_ptr NULL_QUEUE_PTR;

} // namespace kestrelFlow

#endif
