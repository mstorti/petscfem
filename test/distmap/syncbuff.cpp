// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: syncbuff.cpp,v 1.1 2004/01/08 21:10:07 mstorti Exp $
#include <list>
#include <src/distcont.h>

extern int SIZE, MY_RANK;

class TrivialPartitioner {
public:
  int processor(int j) { return int((j*SIZE)/M);};
  void processor(const VT &k,int &nproc,int *plist);
};

void 
TrivialPartitioner::processor(const VT &k,int &nproc,int *plist) {
  nproc=1;
  plist[0] = processor(k.first);
}

template<typename T>
class SyncBuffer : public list<typename T> {
 public:
  
};

