// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: spdirect.h,v 1.4 2002/09/05 18:23:52 mstorti Exp $
#ifndef SPDIRECT_H
#define SPDIRECT_H

#include <petscmat.h>

#include <src/pfmat.h>
#include <src/sparse.h>
#include <src/sparse2.h>

class ThashPFMat : public PFMat {
  TextHashTable thash;
 public:
  void set_option(const char *key,const char *value) {
    thash.add_entry(key,value);
  }
  void get_option(const char *key,const char *&value) const {
    // Remove constness asuming this is OK.
    // (Should fix the `TextHashTable' class
    ((TextHashTable &)thash).get_entry(key,value);
  }
  void get_option(const char *name,int *retval,int defval=0,int n=1) {
    ::get_int(&thash,name,retval,defval,n);
  }
  void get_option(const char *name,double *retval,int defval=0,int n=1) {
    ::get_double(&thash,name,retval,defval,n);
  }
  void get_option(const char *name,string &retval,int defval=0,int
		  n=1 ) {
    get_string(&thash,name,retval,defval,n); }


  /// Not defined yet
  void set_option(const char *name,int *val,int n=1) {
    thash.add_entry(name,val,n); }
  void set_option(const char *name,double *val,int n=1) {
    thash.add_entry(name,val,n); }
};

/// Direct solver. (May be PETSc or SuperLU)
class SparseDirect : public ThashPFMat {
  TextHashTable thash;
  // Does nothing
  int build_sles(TextHashTable *thash,char *name=NULL) {return 0;};
public:
  SparseDirect(int N,char * opt = "PETSc") {
    A_p = Sparse::Mat::dispatch(opt,&thash);
    A_p->resize(N,N);
  }
  /// Returns size
  int size(int j) { return (j==0 ? A_p->rows() : A_p->cols()); }
  /// Not used
  int set_profile_a(int j, int k) { return 0; }
  
  Sparse::Mat *A_p;
  /// destructor
  ~SparseDirect() { A_p->clear(); delete A_p;};
  /// does nothing here (only sequential use...)
  int assembly_begin_a(MatAssemblyType type) { return 0;};
  /// does nothing here (only sequential use...)
  int assembly_end_a(MatAssemblyType type) { return 0;};
  /// Resizes the underlying #A# matrix. ????
  int create_a();
  int set_value_a(int row,int col,Scalar value,
		  InsertMode mode=ADD_VALUES);
  /// Sets all values of the operator to zero.
  int clean_mat_a() {A_p->clear(); return 0;};
  // Does nothing
  int clean_prof_a() {return 0;};
  /** Solve the linear system 
      @param res (input) the rhs vector
      @param dx (input) the solution vector
  */ 
  /// clear memory (almost destructor)
  int clean_factor_a() { return 0;};

  int factor_and_solve_a(Vec &res,Vec &dx);
  /** Solve the linear system 
      @param res (input) the rhs vector
      @param dx (input) the solution vector
  */ 
  int solve_only_a(Vec &res,Vec &dx);
  /// returns the number of iterations spent in the last solve
  int its() {return 1;};
  /// Prints the matrix to a PETSc viewer
  int view(Viewer viewer=NULL) {A_p->print(); return 0;};
  /// Derive this if you want to manage directly the preconditioning. 
  int set_preco(const string & preco_type) {return 0;};
  /// Duplicate matrices (currently not implemented for IISDMat)
  int duplicate_a(MatDuplicateOption op,const PFMat &B);
};

#endif
