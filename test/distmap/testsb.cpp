// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: testsb.cpp,v 1.4 2004/01/18 20:18:35 mstorti Exp $

#include <unistd.h>
#include <list>
#include <iostream>
#include <src/distcont.h>
#include <src/distcont2.h>
#include <src/debug.h>
#include <algorithm>
#include <cassert>
#include <cstdio>

#include <src/syncbuff.h>
#include <src/syncbuff2.h>

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

  int c,m=7,N=10,k=0,ierr=0;
  PetscTruth flg;
  PetscInitialize(&argc,&argv,NULL,NULL);
  MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
  MPI_Comm_rank (MPI_COMM_WORLD, &MY_RANK);

  ierr = PetscOptionsGetInt(PETSC_NULL,"m",&m,&flg);
  CHKERRA(ierr); 

  while ((c = getopt (argc, argv, "m:N:k:")) != -1) {
    switch (c) {
    case 'm':
      sscanf(optarg,"%d",&m);
      break;
    case 'N':
      sscanf(optarg,"%d",&N);
      break;
    case 'k':
      sscanf(optarg,"%d",&k);
      break;
    default:
      PetscPrintf(PETSC_COMM_WORLD,"bad option: \"%c\"\n",c);
      abort();
    }
  }
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] m %d, N %d, k %d\n",MY_RANK,m,N,k);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  PetscFinalize();
  exit(0);
			  
#if 0
  SyncBuffer<PO> sb;

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

  sb.print();

  // sb.check_pack();
#else
  KeyedOutputBuffer kbuff;
  AutoString s;
  
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

#if 1
  FILE *out = fopen("output.dat","w");
  KeyedLine::output = out;
  KeyedLine::print_keys = 0;
  // kbuff.sort_by_key = 0;
  kbuff.flush();
  fclose(out);
#else
  kbuff.flush();
#endif

#endif

  PetscFinalize();
  exit(0);
}
