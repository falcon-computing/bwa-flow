#ifndef QUEUE_H
#define QUEUE_H

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
    bool empty() {
      return data_queue.empty();
    }

    void pop(U &item) {
      while (!data_queue.pop(item)) {
        boost::this_thread::sleep_for(boost::chrono::microseconds(100));
      }
    }

    void push(U item) {
      while (!data_queue.push(item)) {
        boost::this_thread::sleep_for(boost::chrono::microseconds(100));
      }
    }

    bool async_pop(U &item) {
      return data_queue.pop(item);
    }

    bool async_push(U item) {
      return data_queue.push(item);
    }

  private:
    boost::lockfree::queue<U, boost::lockfree::capacity<DEPTH> 
      > data_queue;
};

const boost::shared_ptr<QueueBase> NULL_QUEUE_PTR;

} // namespace kestrelFlow

#endif
