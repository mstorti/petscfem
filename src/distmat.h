// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: distmat.h,v 1.5 2001/08/06 01:07:36 mstorti Exp $
#ifndef DISTMAT_H
#define DISTMAT_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This the row of the matrix. 
*/ 
class Row : public  map<int,double> {
  /// define their own members for packing and unpacking
  int size_of_pack() const;
  /// define their own members for packing and unpacking
  void pack(char *&buff) const;
};

/// This is the basic distributed matrix class. 
typedef DistMap<int,Row> DistMat;

/// This partitioner is based on the dofmap of the mesh. 
class DofmapPartitioner : public Partitioner {
  /// Pointer to the dofmap 
  const Dofmap *dofmap;
public:
  /// Dof partitioning (currently based on ranges of the dof's).  
  int dofpart(int row);
  /// Constructor from the dofmap
  DofmapPartitioner(const Dofmap *dfm);
  /// Destructor
  ~DofmapPartitioner();
};

/** A distributed map<int,int,double> class
 */ 
class DistMatrix : public DistMat {
public:
  /// Specific ``intert'' routine. 
  void insert_val(int i,int j,double v);
  /// Specific function for retrieving values. 
  double val(int i,int j);
  /// Constructor (defines partitioner and communicator). 
  DistMatrix(const Dofmap *dfm, MPI_Comm comm=MPI_COMM_WORLD) 
    : DistMat(new DofmapPartitioner(dfm),comm) {};
  /// Destructor (deletes partitioner).
  ~DistMatrix() {delete part;};
};

#endif
