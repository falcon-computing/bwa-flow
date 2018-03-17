#include <iostream>
#include "TestCommon.h"
#include "MPIChannel.h"

TEST_F(ChannelTests, ChannelSetup) {
  SourceChannel ch1(&link_, 0);
  SourceChannel ch2(&link_, 0);
  SinkChannel   ch3(&link_, 0);

  ASSERT_GT(ch2.getId(), ch1.getId());
  ASSERT_GT(ch3.getId(), ch2.getId());
  ASSERT_GT(ch3.getId(), ch1.getId());

  ASSERT_NE(ch1.getTag(ch1.req), ch2.getTag(ch2.req));
  ASSERT_NE(ch2.getTag(ch2.req), ch3.getTag(ch3.req));
  ASSERT_NE(ch3.getTag(ch3.req), ch1.getTag(ch1.req));
}

static void dispatch_msg(SourceChannel* ch) {
  int rank  = MPI::COMM_WORLD.Get_rank();
  for (int p = 0; p < 10; p++) {
    char msg[5];
    sprintf(msg, "%d", p);
    ch->send(msg, strlen(msg));
  }
  ch->retire();
}

TEST_F(ChannelTests, SourceChannel) {
  // run tests only if there are more than 1 proc
  if (MPI::COMM_WORLD.Get_size() <= 1) {
    return;
  }

  int rank  = MPI::COMM_WORLD.Get_rank();
  int nproc = MPI::COMM_WORLD.Get_size();

  // run tests only if there are more than 3 proc
  ASSERT_GT(nproc, 3);

  SourceChannel ch(&link_, 0);
  boost::thread t;
  if (rank == 0) {
    t = boost::thread(boost::bind(dispatch_msg, &ch));
  }
  while (!ch.recvFinished()) {
    int length;
    char* msg = (char*)ch.recv(length);
    if (length > 0) {
      VLOG(1) << "message " << msg << " for " << rank;
      ASSERT_LT(atoi(msg), 10);
      ASSERT_GE(atoi(msg), 0);
      free(msg);
    }
  }
  if (rank == 0) {
    t.join();
  }
}

static void gather_msg(SinkChannel* ch) {
  // run tests only if there are more than 1 proc
  if (MPI::COMM_WORLD.Get_size() <= 1) {
    return;
  }

  int rank  = MPI::COMM_WORLD.Get_rank();
  int nproc = MPI::COMM_WORLD.Get_size();
  int counter = 0;
  while (!ch->recvFinished()) {
    int length;
    char* msg = (char*)ch->recv(length);

    if (length > 0) {
      VLOG(1) << "message " << msg;

      ASSERT_LT(atoi(msg), 10);
      ASSERT_GE(atoi(msg), 0);

      free(msg);
      counter++;
    }
  }
  ASSERT_EQ(3*nproc, counter);
}

TEST_F(ChannelTests, SinkChannel) {
  // run tests only if there are more than 1 proc
  if (MPI::COMM_WORLD.Get_size() <= 1) {
    return;
  }

  int rank  = MPI::COMM_WORLD.Get_rank();
  int nproc = MPI::COMM_WORLD.Get_size();

  // run tests only if there are more than 3 proc
  ASSERT_GT(nproc, 3);

  SinkChannel ch(&link_, 0);
  boost::thread t;
  if (rank == 0) {
    t = boost::thread(boost::bind(gather_msg, &ch));
  }
  for (int k = 0; k < 3; k++) {
    char msg[10] = {0};
    sprintf(msg, "%d", rank);
    ch.send(msg, strlen(msg));
  }
  ch.retire();
  t.join();
}
