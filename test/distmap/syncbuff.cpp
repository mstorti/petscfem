// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: syncbuff.cpp,v 1.2 2004/01/09 03:08:51 mstorti Exp $
#include <vector>
#include <iostream>
#include <src/distcont.h>
#include <src/distcont2.h>
#include <algorithm>
#include <cassert>

extern int SIZE, MY_RANK;

using namespace std;

class PO  {
public:
  int k;
  PO(int j) : k(j) {}
  void print() { cout << this->k << endl ; }
};

class TrivialPartitioner {
public:
  void processor(const PO &po,int &nproc,int *plist);
};

void 
TrivialPartitioner::processor(const PO &po,int &nproc,int *plist) {
  nproc=1;
  plist[0] = 0;
}

typedef DistCont<vector<PO>,PO,TrivialPartitioner> Container;

// int Container::size_of_pack(const PO &p) const {
//   return sizeof(int)+sizeof(double);
// }

template<typename T>
class SyncBuffer :   
  public DistCont<vector<T>,T,
		  TrivialPartitioner,iter_mode=random_iter_mode> 
  {
    typedef typename DistCont<vector<T>,T,
    TrivialPartitioner,iter_mode=random_iter_mode> cont;
    TrivialPartitioner part;
    public:
    SyncBuffer() : cont(&part) {}
    void sort() { /* not implemented yet */; assert(0); }
    void print() { 
      scatter();
      if (!MY_RANK) {
	typename vector<T>::iterator q;
	for (q=begin(); q!=end(); q++) { q->print(); }
      }
    }
  };

int main() {
  SyncBuffer<PO> sb;
  for (int j=0; j<10; j++) sb.push_back(j);
  sb.print();
}
