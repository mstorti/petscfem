// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: sparse.h,v 1.1 2001/09/20 18:49:40 mstorti Exp $
#ifndef SEQMAT_H
#define SEQMAT_H

#include <map>

#include <fem.h>
#include <utils.h>
#include <dofmap.h>
#include <elemset.h>
#include <pfmat.h>

namespace Sparse {

  typedef vector<int> Indx;
      
  class SeqVec : public map<int,double> {
    int len;
  public:
    SeqVec() {len=0;};
    SeqVec(const SeqVec &v) {*this = v;};
    SeqVec(const Indx &I,const SeqVec &v);
    
    int length();
    //      double & operator[](int j);
    //      const double & operator[](int j) const;
    double get(int j) const;
    void get(const Indx &I,SeqVec &V) const; 

    Seqvec & set(int j,double v);
    Seqvec & set(Indx &I,const SeqVec v);
  }

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
