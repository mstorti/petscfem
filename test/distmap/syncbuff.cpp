// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: syncbuff.cpp,v 1.5 2004/01/09 17:03:16 mstorti Exp $
#include <list>
#include <iostream>
#include <src/distcont.h>
#include <src/distcont2.h>
#include <src/debug.h>
#include <algorithm>
#include <cassert>
#include <cstdio>

int SIZE, MY_RANK;

using namespace std;

const int N=100;
const int m=7;

class PO  {
public:
  int k;
  int *elems;
  int nelem;
  PO() : k(0), elems(NULL), nelem(0) {}
  ~PO() { if(elems) delete[] elems; }
  PO(const PO &po);
  void print();
  friend int operator<(const PO& left, const PO& right);
  int size_of_pack() const;
  void pack(char *&buff) const;
  void unpack(const char *& buff);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int PO::size_of_pack() const { return (2+nelem)*sizeof(int); }

// Copy the int and double to the buffer. Update pointer *buff
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void PO::pack(char *&buff) const {
  memcpy(buff,&k,sizeof(int));
  buff += sizeof(int);
  memcpy(buff,&nelem,sizeof(int));
  buff += sizeof(int);
  memcpy(buff,elems,nelem*sizeof(int));
  buff += nelem*sizeof(int);
}
  
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void PO::unpack(const char *& buff) {
  memcpy(&k,buff,sizeof(int));
  buff += sizeof(int);
  memcpy(&nelem,buff,sizeof(int));
  buff += sizeof(int);
  assert(!elems);
  elems = new int[nelem];
  memcpy(elems,buff,nelem*sizeof(int));
  buff += nelem*sizeof(int);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int operator<(const PO& left, const PO& right) {
  return left.k<right.k;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
PO::PO(const PO &po) {
  if (this==&po) return;
  *this = po;
  elems = new int[nelem];
  for (int j=0; j<nelem; j++) elems[j] = po.elems[j];
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void PO::print() {
  cout << "key: " << k << ", elems: ";
  for (int j=0; j<nelem; j++) cout << elems[j] << " ";
  cout << endl;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class TrivialPartitioner {
public:
  void processor(const PO &po,int &nproc,int *plist);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void 
TrivialPartitioner::processor(const PO &po,int &nproc,int *plist) {
  nproc=1;
  plist[0] = 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
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
  void check_pack();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int DistCont<list<PO>,PO,TrivialPartitioner>
::size_of_pack(const PO &po) const {
  return po.size_of_pack();
}

// Copy the int and double to the buffer. Update pointer *buff
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistCont<list<PO>,PO,TrivialPartitioner>
::pack(const PO &po,char *&buff) const { po.pack(buff); }
  
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistCont<list<PO>,PO,TrivialPartitioner>
::unpack(PO &po,const char *& buff) { po.unpack(buff); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistCont<list<PO>,PO,TrivialPartitioner>
::combine(const PO &po) { 
  push_back(po);
#if 0
  back() = po;
  back().elems = new int[po.nelem];
  for (int j=0; j<po.nelem; j++) back().elems[j] = po.elems[j];
#endif
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<typename T>
void SyncBuffer<T>::check_pack() {
  print();
  int bufsize = 0;
  typename list<T>::iterator q;
  for (q=begin(); q!=end(); q++) 
    bufsize += q->size_of_pack();
  char *buff = new char[bufsize], *p=buff, *pe = buff+bufsize;
  printf("%d elems in buffer\n",list<T>::size());
  for (q =begin(); q!=end(); q++) q->pack(p);
  assert(p == pe);
  clear();
  const char *pp = buff;
  int nelems=0;
  while (pp<pe) {
    nelems++;
    push_back();
    back().unpack(pp);
  }
  assert(pp==pe);
  print();
  printf("%d elems recovered from buffer\n",list<T>::size());
}

int main(int argc,char **argv) {

  PetscInitialize(&argc,&argv,NULL,NULL);
  MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
  MPI_Comm_rank (MPI_COMM_WORLD, &MY_RANK);

  Debug debug;
  debug.init();
  //  debug.activate();
  SyncBuffer<PO> sb;
  int N=10;
  for (int j=0; j<N; j++) {
    sb.push_back();
    int k = SIZE*j+MY_RANK;
    sb.back().k = k;
    int roof = k - (k % m) +m;
    int nelem = roof-k;
    sb.back().nelem = nelem;
    int *elems = new int[nelem];
    sb.back().elems = elems;
    for (int jj=0; jj<nelem; jj++) elems[jj] = k+jj;
    sb.back().print();
  }

  // sb.check_pack();

  debug.trace("antes de sb.print()");
  sb.print();
  PetscFinalize();
  exit(0);
}
