// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: distmat.h,v 1.4 2001/08/03 17:07:25 mstorti Exp $
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
  DistMatrix(const Dofmap *dfm, MPI_Comm comm=MPI_COMM_WORLD) 
    : DistMat(new DofmapPartitioner(dfm),comm) {};
  ~DistMatrix() {delete part;};
};

#endif
