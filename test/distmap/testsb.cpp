// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: testsb.cpp,v 1.1 2004/01/11 16:22:09 mstorti Exp $
#include <list>
#include <iostream>
#include <src/distcont.h>
#include <src/distcont2.h>
#include <src/debug.h>
#include <algorithm>
#include <cassert>
#include <cstdio>

#include <src/syncbuff.h>

int SIZE, MY_RANK;

using namespace std;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
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

SYNC_BUFFER_FUNCTIONS(PO);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int PO::size_of_pack() const { return (2+nelem)*sizeof(int); }

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
  for (int j=0; j<nelem; j++) cout << elems[j] << " ";
  cout << endl;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int main(int argc,char **argv) {

  PetscInitialize(&argc,&argv,NULL,NULL);
  MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
  MPI_Comm_rank (MPI_COMM_WORLD, &MY_RANK);

  Debug debug;
  debug.init();
  //  debug.activate();
#if 0
  SyncBuffer<PO> sb;
  int N=10;
  int m=7;

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

  debug.trace("antes de sb.print()");
  sb.print();

  // sb.check_pack();
#else
  KeyedOutputBuffer kbuff;
  AutoString s;
  
  int N=10;
  int m=7;
  for (int j=0; j<N; j++) {
    int k = SIZE*j+MY_RANK;

    int roof = k - (k % m) +m;
    int nelem = roof-k;

    s.clear();
    for (int jj=0; jj<nelem; jj++) s.cat_sprintf("%d ",k+jj);
    kbuff.push(k,s);
    // kbuff.back().print();
  }
  // kbuff.check_pack();

#if 0
  FILE *out = fopen("output.dat","w");
  KeyedLine::output = out;
  KeyedLine::print_line_numbers = 0;
  kbuff.print();
  fclose(out);
#else
  kbuff.print();
#endif

#endif

  PetscFinalize();
  exit(0);
}
