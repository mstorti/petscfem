// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: idmap.h,v 1.3 2001/05/30 03:58:50 mstorti Exp $
 
#ifndef IDMAP_H
#define IDMAP_H

#include <map>
#include <set>

#include "fastlib.h"

enum map_type {NULL_MAP, IDENTITY_MAP};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** @name idmap package*/
//@{
/// Stores rows of the Q matrix 
typedef map<int,double> row_t;

struct IdMapEntry {
  int j;
  double coef;
};

typedef FastVector<IdMapEntry> IdMapRow;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** `daxpy' type operation ( y <- y + a*x ) for `idmap' rows. 
    @author M. Storti
    @param y (input/output ) row to be modified. 
    @param a (input) scalar coeficient. 
    @param x (input) row to be added. 
*/ 
void axpy(row_t &y,double const a,const row_t &x);

void print(const row_t &y);

/// Stores columns of the Q matrix 
typedef set<int> col_t;

typedef map<int,row_t*> row_map_t;
typedef map<int,col_t*> col_map_t;


//@}
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/**
``idmap'' class stores sparse m x n matrices that are ``close'' to a
permutation. 
*/
class idmap {
public:

  ~idmap();

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Constructor from an initial permutation given by ident and iident.
      @author M. Storti
      @param m row dimension of Q
      @param n column dimension of Q
      @param ident integer array (permutation). 
      @param iident integer array (permutation). 
   */ 
  idmap(int m,int n,int *ident,int *iident);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Constructor for a square (m=n) map to identity or null.
      @author M. Storti
      @param m (input) dimension of matrix
      @param ival (input) may be IDENTITY\_MAP or NULL\_MAP
   */ 
  idmap(int m,map_type mt);
  /// 

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Sets element i,j to val.
      @author M. Storti
      @param i row index
      @param j column index
      @param val coefficient to be entered at position i,j. (default=0)
   */ 
  void set_elem(const int i,const int j,const double val=0.);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Deletes  the i-th row.
      @author M. Storti
      @param i index of row to be deleted
   */ 
  void del_row(const int i);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Sets the i-th row to row.
      @author M. Storti
      @param i row index
      @param row row to be entered in the i-th position
   */ 
  void row_set(const int i,const row_t &row);

  /** Gets the i-th row.
      @author M. Storti
      @param i index of row to be retrieved
      @param row retrieved row (output)
   */ 
  void get_row(const int i,row_t &row);

  /// overloaded in order to run faster
  void get_row(const int i,IdMapRow &row);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Gets the j-th column to col. 
      @author M. Storti
      @param j column index
      @param col column to be insert in j-th position
   */ 
  void get_col(const int j,row_t &col);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Gets the i,j value Q\_{ij}.
      @author M. Storti
      @param i row index
      @param j column index
      @param val coefficient Q\_{ij}
  */ 
  void get_val(const int i,const int j,double &val);
  
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Deletes a column. 
      @author M. Storti
      @param j (input) column indx
  */ 
  void del_col(const int j);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Sets the j-th column to col.
      @author M. Storti
      @param j (input) column index
      @param col (input) column to be entered in the j-th position
   */ 
  void column_set(const int j,row_t &col);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Prints the Q matrix by rows (in sparse form).
      @author M. Storti
      @param s optional string
   */
  void print(char *s = NULL);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Prints the Q matrix by cols. 
      @author M. Storti
   */

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  void print_by_col(char *s=NULL);
  /** Checks internal consistency of idmap. 
      Checks that any entry (by row) has the corresponding entry by
      column, and viceversa. Also, that all entries in `special rows'
      are non-zero. 
      @author M. Storti 
      @param s optional string
   */
  void check();

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Solves an equation of the form Q x = y for x. 
      Assumes that rank(Q) = n = dim(x). 
      @author M. Storti 
  */
  void solve(double *x,double *y);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /**Extracts the matrix connected to a given column `j'. 
     This routine is the key part for solving a linear system
     associated with the underlying matrix. (In a future I will write a
     similar one for index  row).
     @author M. Storti 
     @param j (input) column index to get the connected matrix
     @param iindx set of row indices connected to j
     @param jindx set of column indices  connected to j
     @param q Newmat matrix returned. Entry (k,l) in q
     corresponds to the k-th entry in iindx and the l-th entry in jindx
  */
  void get_block_matrix(int j,set<int> &iindx,set<int> &jindx,Matrix &qq);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Redefines the matrix through a reordering of the columns. 
      @author M. Storti
      @param perm (input) The vector that remaps columns. 
      new\_j\_indx = perm[old\_j\_indx-1]
  */ 
  void remap_cols(int *perm);
  void remap_cols(int *perm,idmap &idnew);
  
  void print_statistics(void) const;

private:
  /// row and column dimension
  int m,n;
  /// pointers row <-> column
  int *ident,*iident;
  /// for `special' rows, maps row index to a pointer to a `special row'
  row_map_t *row_map;
  /// for `special' columns, maps column index to a pointer `special column'
  col_map_t *col_map;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** @name Class `idmap' utility functions */
//@{
/// prints a column with an optional string
void print(col_t &col,char *s=NULL);

/// prints a row with an optional string
void print(row_t &row,char *s=NULL);

/// voids a row
void erase_null(row_t &row);
//@}

#endif
