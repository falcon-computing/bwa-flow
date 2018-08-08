#ifndef TASK_H
#define TASK_H

class Task {
 public:
  virtual void start(Task* prev_task) {;}
  virtual void finish() = 0;

};
#endif

