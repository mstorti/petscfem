// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: syncbuff.cpp,v 1.4 2004/01/09 15:00:42 mstorti Exp $
#include <list>
#include <iostream>
#include <src/distcont.h>
#include <src/distcont2.h>
#include <src/debug.h>
#include <algorithm>
#include <cassert>

int SIZE, MY_RANK;

using namespace std;

class PO  {
public:
  int k;
  PO(int j=0) : k(j) {}
  void print() { cout << this->k << endl ; }
  int less(PO p2) { return this->k<p2.k; }
  friend int operator<(const PO& left, const PO& right);
};

int operator<(const PO& left, const PO& right) {
  return left.k<right.k;
}

class TrivialPartitioner {
public:
  void processor(const PO &po,int &nproc,int *plist);
};

void 
TrivialPartitioner::processor(const PO &po,int &nproc,int *plist) {
  nproc=1;
  plist[0] = 0;
}

template<typename T>
class SyncBuffer :   
  public DistCont<list<T>,T,TrivialPartitioner>   {
  //    typedef typename DistCont<list<T>,T,TrivialPartitioner> cont;
  TrivialPartitioner part;
public:
  SyncBuffer() : 
    DistCont<list<T>,T,TrivialPartitioner>(&part) {}
  void print() { 
    scatter();
    sort();
    if (!MY_RANK) {
      typename list<T>::iterator q;
      for (q=begin(); q!=end(); q++) { q->print(); }
    }
  }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int DistCont<list<PO>,PO,TrivialPartitioner>
::size_of_pack(const PO &po) const {
  return sizeof(PO);
}

// Copy the int and double to the buffer. Update pointer *buff
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistCont<list<PO>,PO,TrivialPartitioner>
::pack(const PO &po,char *&buff) const {
  memcpy(buff,&po.k,sizeof(int));
  buff += sizeof(int);
}
  
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistCont<list<PO>,PO,TrivialPartitioner>
::unpack(PO &po,const char *& buff) {
  memcpy(&po.k,buff,sizeof(int));
  buff += sizeof(int);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistCont<list<PO>,PO,TrivialPartitioner>
::combine(const PO &po) { push_back(po); };
 
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool PO_less(PO po1,PO po2) { return po1.k<po2.k; }

int main(int argc,char **argv) {

  PetscInitialize(&argc,&argv,NULL,NULL);
  MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
  MPI_Comm_rank (MPI_COMM_WORLD, &MY_RANK);

  Debug debug;
  debug.init();
  debug.activate();
  SyncBuffer<PO> sb;
  int N=10;
  for (int j=0; j<N; j++) sb.push_back(SIZE*j+MY_RANK);
  debug.trace("antes de sb.print()");
  sb.print();
  PetscFinalize();
  exit(0);

}
