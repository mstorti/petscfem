// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: distmat.h,v 1.3 2001/08/02 19:50:22 mstorti Exp $
#ifndef DISTMAT_H
#define DISTMAT_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Test for the distributed map class
// A distributed map<int,int,double> class
class Row : public  map<int,double> {
  int size_of_pack() const;
  void pack(char *&buff) const;
};

typedef DistMap<int,Row> DistMat;
#if 0
class DistMat : public DistMap<int,Row> {
public:
  DistMat(Partitioner *p) : DistMap<int,Row>(p) {};
};
#endif
//typedef map<int,Row> BasMap;
//DistMat::DistMat(const Dofmap *dfm) : DistMap<int,Row>(new DofmapPartitioner(dfm);) {};

class DofmapPartitioner : public Partitioner {
  const Dofmap *dofmap;
public:
  int dofpart(int row);
  DofmapPartitioner(const Dofmap *dfm);
  ~DofmapPartitioner();
};

class DistMatrix : public DistMat {
public:
  void insert_val(int i,int j,double v);
  double val(int i,int j);
  DistMatrix(const Dofmap *dfm) : DistMat(new DofmapPartitioner(dfm)) {};
  ~DistMatrix() {delete part;};
};

#endif
