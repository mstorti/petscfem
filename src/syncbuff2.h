// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: syncbuff2.h,v 1.4 2005/10/19 17:40:33 mstorti Exp $
#ifndef PETSCFEM_SYNCBUFF2_H
#define PETSCFEM_SYNCBUFF2_H

#include <list>
#include <iostream>
#include <src/distcont.h>
#include <src/distcont2.h>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <src/autostr.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<typename T>
void TrivialPartitioner<T>
::processor(const T &t,int &nproc,int *plist) {
  nproc=1;
  plist[0] = 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<typename T>
void SyncBuffer<T>::print() { 
  // Does the scatter
  this->scatter(); 
  // Does the sorting
  if (sort_by_key) this->sort();
  // Prints all elements in master
  if (!MY_RANK) {
    typename list<T>::iterator q;
    for (q = this->begin(); q != this->end(); q++) { q->print(); }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<typename T>
void SyncBuffer<T>::check_pack() {
  // This is very useful to debug, since it's very comon to make
  // errors in the pack/unpack routines.  As the sequential version
  // doesn't need to do the scatter one has to debug the packing and
  // unpacking routines onlw in parallel.  This routine simulates
  // packing and unpacking, checking that total sizes are respected.
  // To use, generate a container with elements and apply.  It prints
  // the elements, packs all of them in a buffer, cleans the container
  // and unpacks all the elements from it.
  print();

  // Computes the size of the buffer by calling #size_of_pack# 
  // for all elements. 
  int bufsize = 0;
  typename list<T>::iterator q;
  for (q = this->begin(); q != this->end(); q++) 
    bufsize += q->size_of_pack();

  // Create the buffer
  char *buff = new char[bufsize], *p=buff, *pe = buff+bufsize;
  printf("%d elems in buffer\n",list<T>::size());

  // Pack all elements in buffer
  for (q = this->begin(); q != this->end(); q++) q->pack(p);
  assert(p == pe);

  // Clean the container
  this->clear();

  // Unpacks all elements from the buffer into
  // the container. 
  const char *pp = buff;
  int nelems=0;
  while (pp<pe) {
    nelems++;
    this->push_back();
    this->back().unpack(pp);
  }
  assert(pp==pe);
  print();

  printf("%d elems recovered from buffer\n",list<T>::size());
}

#endif
