// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: distmat.h,v 1.1 2001/08/01 22:55:33 mstorti Exp $
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
typedef map<int,Row> BasMap;

class DistMatrix : public DistMat {
public:
  void insert_val(int i,int j,double v);
  double val(int i,int j);
};

#endif
