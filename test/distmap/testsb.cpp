// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: testsb.cpp,v 1.8 2005/01/26 10:54:48 mstorti Exp $

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
class KeyedObject  {
public:
  int k;
  int *elems;
  int nelem;
  KeyedObject() : k(0), elems(NULL), nelem(0) {}
  ~KeyedObject() { if(elems) delete[] elems; }
  KeyedObject(const KeyedObject &po);
  void print();
  friend int operator<(const KeyedObject& left, const KeyedObject& right);
  int size_of_pack() const;
  void pack(char *&buff) const;
  void unpack(const char *& buff);
};

SYNC_BUFFER_FUNCTIONS(KeyedObject);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int KeyedObject::size_of_pack() const { return (2+nelem)*sizeof(int); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedObject::pack(char *&buff) const {
  memcpy(buff,&k,sizeof(int));
  buff += sizeof(int);
  memcpy(buff,&nelem,sizeof(int));
  buff += sizeof(int);
  memcpy(buff,elems,nelem*sizeof(int));
  buff += nelem*sizeof(int);
}
  
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedObject::unpack(const char *& buff) {
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
int operator<(const KeyedObject& left, const KeyedObject& right) {
  return left.k<right.k;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
KeyedObject::KeyedObject(const KeyedObject &po) {
  if (this==&po) return;
  *this = po;
  elems = new int[nelem];
  for (int j=0; j<nelem; j++) elems[j] = po.elems[j];
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedObject::print() {
  cout << "k: " << k << ", ";
  for (int j=0; j<nelem; j++) 
    cout << elems[j] << " ";
  cout << endl;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int main(int argc,char **argv) {

  // Checks the Syncbuffer (k=0) and KeyedOutputBuffer (k=1) classes.
  // N rows of integers are loaded. Row `j' goes from `j' to
  // `k - (k % m) +m' i.e. to the next integer multiple of `m'.
  // For instance, for `m=7, N=15' the following rows are loaded

  // row 0: 0 1 2 3 4 5 6 
  // row 1: 1 2 3 4 5 6 
  // row 2: 2 3 4 5 6 
  // row 3: 3 4 5 6 
  // row 4: 4 5 6 
  // row 5: 5 6 
  // row 6: 6 
  // row 7: 7 8 9 10 11 12 13 
  // row 8: 8 9 10 11 12 13 
  // row 9: 9 10 11 12 13 
  // row 10: 10 11 12 13 
  // row 11: 11 12 13 
  // row 12: 12 13 
  // row 13: 13 
  // row 14: 14 15 16 17 18 19 20 

  // Row `j' is loaded in processor `j % SIZE', where `SIZE' is the
  // number of processors. Then the rows are scattered to the master
  // processor and printed.

  // For `k=0' a `KeyedObject' class is defined which is a dynamic
  // integer vector with a print function. For the `KeyedOutputBuffer'
  // class the integers are printed in the text line with the help of
  // an `AutoString' object. 

  // Arguments:
  //  N: number of rows
  //  m: rows go until the next multiple of `m'
  //  k: check SyncBuffer class (k=0) or
  //             `KeyedOutputBuffer' class (k=1)
  //  s: sort records by key (only for `k=1')
  //  p: print keys in `KeyedOutputBuffer'

  int k,output=0,m=7,N=10,keyed=0,ierr=0,
    sort_by_key=1,print_keys=1,print_newlines=1;

  PetscTruth flg;
  PetscInitialize(&argc,&argv,NULL,NULL);
  MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
  MPI_Comm_rank (MPI_COMM_WORLD, &MY_RANK);

  ierr = PetscOptionsGetInt(PETSC_NULL,"-m",&m,&flg);
  CHKERRA(ierr); 
  ierr = PetscOptionsGetInt(PETSC_NULL,"-N",&N,&flg);
  CHKERRA(ierr); 
  ierr = PetscOptionsGetInt(PETSC_NULL,"-k",&keyed,&flg);
  CHKERRA(ierr); 
  ierr = PetscOptionsGetInt(PETSC_NULL,"-o",&output,&flg);
  CHKERRA(ierr); 
  ierr = PetscOptionsGetInt(PETSC_NULL,"-s",&sort_by_key,&flg);
  CHKERRA(ierr); 
  ierr = PetscOptionsGetInt(PETSC_NULL,"-p",&print_keys,&flg);
  CHKERRA(ierr); 
  ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&print_newlines,&flg);
  CHKERRA(ierr); 

  if (!keyed) {
    PetscPrintf(PETSC_COMM_WORLD,
		"Testing basic SyncBuffer class\n");
		// "Disordered elements: \n");
    SyncBuffer<KeyedObject> sb;
    KeyedObject t;
    
    k=MY_RANK;
    while (k<N) {
      sb.push_back(t);
      sb.back().k = k;
      int roof = k - (k % m) +m;
      int nelem = roof-k;
      sb.back().nelem = nelem;
      int *elems = new int[nelem];
      sb.back().elems = elems;
      for (int jj=0; jj<nelem; jj++) elems[jj] = k+jj;
      // sb.back().print();
      k += SIZE;
    }

    sb.print();

  } else {

    PetscPrintf(PETSC_COMM_WORLD,"Testing KeyedOutputBuffer class\n");

    KeyedOutputBuffer kbuff;
  
    k=MY_RANK;
    while (k<N) {
      int roof = k - (k % m) +m;
      int nelem = roof-k;

      for (int jj=0; jj<nelem; jj++) kbuff.cat_printf("%d ",k+jj);
      if (!print_newlines) kbuff.cat_printf(" (explicit new-line)\n");
      kbuff.push(k);
      k += SIZE;
    }

    kbuff.sort_by_key = sort_by_key;
    KeyedLine::print_keys = print_keys;
    KeyedLine::print_newlines = print_newlines;
    if (output) {
      FILE *out = fopen("testsb.output.tmp","w");
      KeyedLine::output = out;
      kbuff.flush();
      fclose(out);
    } else kbuff.flush();

  }
  PetscFinalize();
  exit(0);
}
