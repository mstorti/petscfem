// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: sparse.h,v 1.2 2001/09/20 20:13:24 mstorti Exp $
#ifndef SEQMAT_H
#define SEQMAT_H

#include <cstdio>

#include <map>
#include <vector>

namespace Sparse {

  typedef vector<int> Indx;
  typedef map<int,double>::iterator VecIt;
  typedef map<int,double>::const_iterator VecCIt;
  typedef pair<int,double> VecP;
      
  /// A simple sparse vector class (0 indexed). 
  class Vec : public map<int,double> {
    /// Length of the sparse vector
    int len;
    /// Flag indicating where you can add values past the specified length or not 
    int grow_m; 
  public:
    
    /// Constructor from the length
    Vec(int l=0) : grow_m(1) {len=l;};
    /// Return length of the vector 
    int length() {return len;};
    /// Constructor from another vector
    Vec(const Vec &v) {*this = v;};
    /// Insert contents of vector v at position I
    Vec(const Indx &I,const Vec &v);
    /// Get element at specified position
    double get(int j) const;
    /// Get a subvector of elements at position I
    void get(const Indx &I,Vec &V) const; 
    /// Set element at position j
    Vec & set(int j,double v);
    /// Set elements at subvector at position I
    Vec & set(Indx &I,const Vec v);
    /// print elements
    void print();
    /// Set mode if can grow automatically or not
    Vec & grow(int g) { grow_m=g; return *this;};

  };

#if 0
    // Simple sparse matrix class. 
  class SeqMat : public PFMat, public map< int, Row >  {
    typedef map<int,Row>::iterator RowIt;
    typedef map<int,double>::iterator ColIt;
    typedef pair<int,Row> row_p;
    typedef pair<int,double> col_p;

    int nrows,ncols;
    /// PETSc error code 
    int ierr;
  public:
    /// Destructor (calls almost destructor)
    ~SeqMat() {clear();};
    /// clear memory (almost destructor)
    // void clear() {clear();};
    /// Constructor 
    SeqMat() : PFMat() {};
    /** Call the assembly begin function for the underlying PETSc matrices
	@param type (input) PETSc assembly type
	@return A PETSc error code 
    */ 
    int assembly_begin(MatAssemblyType type) { return 0;};
    /** Call the assembly end function for the underlying PETSc matrices
	@param type (input) PETSc assembly type
	@return A PETSc error code 
    */ 
    int assembly_end(MatAssemblyType type) { return 0;};
    /** Sets individual values on the operator #A(row,col) = value#
	@param row (input) first index
	@param col (input) second index
	@param value (input) the value to be set
	@param mode (input) either #ADD_VALUES# (default) or #INSERT_VALUES#
    */ 
    void set_value(int row,int col,Scalar value,InsertMode mode=ADD_VALUES);
    /// Sets all values of the operator to zero.
    int zero_entries() {clear();};
    /** Does nothing for this class.
	Creates the matrix from the profile computed in #da#
	@param da (input) dynamic array containing the adjacency matrix
	of the operator
	@param dofmap (input) the dofmap of the operator (contains
	information about range of dofs per processor. 
	@param debug_compute_prof (input) flag for debugging the process
	of building the operator.
    */ 
    void create(Darray *da,const Dofmap *dofmap_,int debug_compute_prof=0);
    /// Duplicate matrices 
    int duplicate(MatDuplicateOption op,PFMat &A);
    int view(Viewer viewer);
  };
#endif

}

#endif
