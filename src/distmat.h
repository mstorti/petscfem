// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: distmat.h,v 1.6 2001/08/13 00:12:38 mstorti Exp $
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

/// This partitioner is based on the dofmap of the mesh. 
class DofmapPartitioner : public Partitioner {
  /// Pointer to the dofmap 
  const Dofmap *dofmap;
public:
  /// Dof partitioning (currently based on ranges of the dof's).  
  int processor(int row);
  /// Constructor from the dofmap
  DofmapPartitioner(const Dofmap *dfm);
  /// Destructor
  ~DofmapPartitioner();
};

/// This is the basic distributed matrix class. 
typedef DistMap<int,Row,DofmapPartitioner> DistMat;

/** A distributed map<int,int,double> class
 */ 
class DistMatrix : public DistMat {
public:
  /// Specific ``insert'' routine. 
  void insert_val(int i,int j,double v);
  /// Specific function for retrieving values. 
  double val(int i,int j);
  /// Constructor (defines partitioner and communicator). 
  DistMatrix(MPI_Comm comm=MPI_COMM_WORLD) 
    : DistMat(p,comm) {};
  /// Destructor (deletes partitioner).
  ~DistMatrix() {delete part;};
};

#endif
