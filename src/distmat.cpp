/*__INSERT_LICENSE__*/
// $Id: distmat.cpp,v 1.6 2001/08/06 01:07:36 mstorti Exp $
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <petsc.h>

#include <utils.h>
#include <distmap.h>
#include <buffpack.h>
#include <maximizr.h>
#include <distmat.h>

// eckel:=  
// This is required!! See `Thinking in C++, 2nd ed. Volume 1', 
// by Bruce Eckel (http://www.MindView.net)
// Chapter 15: `Polymorphism &  Virtual Functions', 
// paragraph `Pure virtual destructors' 
Partitioner::~Partitioner() {};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int DistMat::
size_of_pack(const DistMat::iterator iter) const {
  int n = iter->second.size();
  // size + row number + size*(int+double)
  return (n+2)*sizeof(int)+n*sizeof(double);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistMat::
pack(const int &k,const Row &row,char *&buff) const {
  int n = row.size();
  Row::const_iterator iter;
  // number of elements in the row
  BufferPack::pack(n,buff);
  // the number of the row (key)
  BufferPack::pack(k,buff);
  // pack the row, each pair (int,double) at a time
  for (iter=row.begin(); iter!=row.end(); iter++) {
    BufferPack::pack(iter->first,buff);
    BufferPack::pack(iter->second,buff);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void DistMat::unpack(int &,Row &,const char *&)"
void DistMat::
unpack(int &k,Row &row,const char *&buff) {
  int n,j,key;
  double val;
  // Clear the row
  row.clear();
  // unpack the number of elements
  BufferPack::unpack(n,buff);
  // unpack the key (number of row)
  BufferPack::unpack(k,buff);
  // unpack the pair (key vals) each at a time
  for (j=0; j<n; j++) {
    BufferPack::unpack(key,buff);
    BufferPack::unpack(val,buff);
    // insert the pair key,val in the row.
    row[key] += val;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int DistMat::processor(const DistMatrix::iterator ) const"
int DistMat::processor(const DistMat::iterator k) const {
  // The number of processor is computed in the `Partitioner' 
  return part->dofpart(k->first);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistMat::
combine(const pair<int,Row> &p) {
  // Insert the pair (int,row) in the matrix
  int n,j;
  map<int,Row>::iterator iter = find(p.first);
  Row::iterator r;
  Row::const_iterator q;
  if (iter == end()) {
    // insert(p); // this doesn't work
    // If entry doesn't exist for this key, then insert a new row for
    // this key 
    (*this)[p.first] = p.second;
  } else {
    // merge the elements of this row to the existing ones
    // existing row
    Row &oldr = iter->second;
    // new row
    const Row &newr = p.second;
    // new element in the row
    for (q = newr.begin(); q!=newr.end(); q++) {
      // r iterator to that key in the old row
      r = oldr.find(q->first);
      if (r == oldr.end()) {
	// the col is not in the row
	// insert a new entry
	oldr.insert(*q);
      } else {
	// add to the existing 
	r->second += q->second;
      }
    }
  }
}

void DistMatrix::insert_val(int i,int j,double v) {
  DistMat::iterator I = find(i);
  Row::iterator J;
  if (I == end()) {
    insert(pair<int,Row>(i,Row()));
    I = find(i);
  }
  Row &row = I->second;
  J = row.find(j);
  if (J == row.end()) {
    row.insert(pair<int,double>(j,v));
  } else {
    J->second += v;
  }
}

double DistMatrix::val(int i,int j) {
  Row::iterator J;
  DistMat::iterator I = find(i);
  if (I == end()) return 0.;
  Row &row = I->second;
  J = row.find(j);
  if (J == row.end()) {
    return 0.;
  } else {
    return J->second;
  }
}  

DofmapPartitioner::
DofmapPartitioner(const Dofmap *dfm) :  dofmap(dfm) {};

int DofmapPartitioner::dofpart(int row) {
  const int *startproc = dofmap->startproc;
  const int *neqproc = dofmap->neqproc;
  int proc;
  for (proc = 0; proc < dofmap->size; proc++) 
    if (row < startproc[proc]+neqproc[proc]) 
      return proc;
}

DofmapPartitioner::~DofmapPartitioner() {};
