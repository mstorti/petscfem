/*__INSERT_LICENSE__*/
// $Id: distmat.cpp,v 1.3 2001/08/02 19:50:22 mstorti Exp $
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <petsc.h>

#include <utils.h>
#include <distmap.h>
#include <buffpack.h>
#include <maximizr.h>
#include <distmat.h>

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
  BufferPack::pack(n,buff);
  BufferPack::pack(k,buff);
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
  BufferPack::unpack(n,buff);
  BufferPack::unpack(k,buff);
  for (j=0; j<n; j++) {
    BufferPack::unpack(key,buff);
    BufferPack::unpack(val,buff);
    row[key] += val;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int DistMat::processor(const DistMatrix::iterator ) const"
int DistMat::processor(const DistMat::iterator k) const {
  return part->dofpart(k->first);
}
#if 0
int DistMat::processor(const DistMat::iterator k) const {
  const int &row = k->first;
  const int *startproc = dofmap->startproc;
  const int *neqproc = dofmap->neqproc;
  int proc;
  for (proc = 0; proc < size; proc++) 
    if (row < startproc[proc]+neqproc[proc]) 
      return proc;
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistMat::
combine(const pair<int,Row> &p) {
  int n,j;
  map<int,Row>::iterator iter = find(p.first);
  Row::iterator r;
  Row::const_iterator q;
  if (iter == end()) {
    printf("[%d] inserting row %d\n",myrank,p.first);
    // insert(p);
    (*this)[p.first] = p.second;
  } else {
    Row &oldr = iter->second;
    const Row &newr = p.second;
    for (q = newr.begin(); q!=newr.end(); q++) {
      r = oldr.find(q->first);
      if (r == oldr.end()) {
	r->second += q->second;
      } else {
	oldr.insert(*q);
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
