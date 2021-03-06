// -*- mode: c++ -*-
/*__INSERT_LICENSE__*/
//$Id merge-with-petsc-233-55-g52bd457 Fri Oct 26 13:57:07 2007 -0300$
#ifndef UTIL2_H
#define UTIL2_H

#include <vector>
#include <newmat.h>
#include <petscksp.h>
#include <src/fstack.h>
// Where are they defined??
#undef min
#undef max
#include <src/fastmat2.h>

typedef double scalarfun(double,void *);

typedef ColumnVector vectorfun(ColumnVector,void *);

int non_symm_eigenvals(const Matrix &A,Matrix &lambda,Matrix &Vre,Matrix &Vim);

int mat_function(const Matrix &A,Matrix &funA,scalarfun fun,void *);

int mat_function(const Matrix &A,Matrix &funA,vectorfun fun,void *);

int vector_divide(Vec &res,Vec a_mass);

void show_mallinfo (char *s=NULL);

#define TRC {PFEM_TRACE(""); show_mallinfo();}
#define TRCS(s) {PFEM_TRACE(s); show_mallinfo();}

#include <time.h>
#include <sys/time.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Manage chronometers.
    @author M. Storti
*/ 
/// 
class Chrono {
private:
  /// structures for `libc' time function calls.
  clock_t start_time;
public:
#if 1
  /// create a chronometer
  Chrono() : start_time(clock()) { }
  /// destroy chronometer
  ~Chrono() {};
#endif
  /// return elapsed CPU time from start
  double elapsed() const {return ((double) (clock() - start_time)) / CLOCKS_PER_SEC;};
  /// reset start time to actual time
  void start() {start_time = clock();};
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** High precision cronometer
    @author M. Storti
*/ 
/// 
class HPChrono {
  double gettod() const {
    struct timeval tv;
    gettimeofday(&tv,0);
    return tv.tv_sec + 1e-6 * tv.tv_usec;
  }
private:
  /// structures for `libc' time function calls.
  double start_time;
public:
  HPChrono() : start_time(0) { }
  /// return elapsed CPU time from start
  double elapsed() const {return gettod()-start_time;};
  /// reset start time to actual time
  void start() {start_time = gettod();};
};

class TextHashTable;
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Reads a hash table from a filestack. 
    @author M. Storti
    @param fstack (input) the file stack from where the hash will be read
    @param thash  (output) the text hash table that is read
*/ 
int read_hash_table(FileStack *& fstack,TextHashTable *& thash);

double pw4(double x);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Reads a series of doubles from a string into a
    vector<double>. 
    Useful when reading hash-tables with a number of double
    arguments. 
    @author M. Storti
    @param v (output) the array of doubles returned
    @return s (input) the string where the doubles are read. 
*/ 
void read_double_array(vector<double> &v,const char * s);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Reads a series of ints from a string into a
    vector<int>. 
    Useful when reading hash-tables with a number of int
    arguments. 
    @author M. Storti
    @param v (output) the array of ints returned
    @return s (input) the string where the ints are read. 
*/ 
void read_int_array(vector<int> &v,const char * s);

/** Equivalent to 'pow' but for integer powers.
    (May be more efficient).
    @param (input)
    @return a reference to the matrix.
*/ 
double int_pow(double base,int exp);

/** Cyclic `rem'. #crem(n,m)# computes the remainder of division 
    of `n' by `m' but (in contrast with the C builtin % operator) 
    is extended cyclically to negative numbers so that for instance, for
    the numbers #n# -4 to 4 we have: #rem (*,3) =
    {-1,0,-2,-1,0,1,2,0,1}# whereas #crem (*,3) = {2,0,1,2,0,1,2,0}#
    [Note: this is duplicated with `modulo' defined in `utils.h' ]
*/
int crem(int j, int m);

/** @name Wrapper to PETSc destroy functions 
    @doc Check if the pointer is NULL.
*/
//@{
/// VecDestroy wrapper
int VecDestroy_maybe(Vec &v);
/// MatDestroy wrapper
int MatDestroy_maybe(Mat &v);
/// KSPDestroy wrapper
int KSPDestroy_maybe(KSP &v);
//@}

#define DELETE_SCLR(name)			\
  if (name) { delete name; name=NULL; }

#define DELETE_VCTR(name)			\
  if (name) { delete[] name; name=NULL; }

#define DELETE_FUNC(delfunc, name)		\
  if (name) { delfunc(name); name=NULL; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class MakeTangentSpace {
private:
  int ndim;
  FastMat2 z;
public:
  void init(int ndim_a);
  void make_tangent(const FastMat2 &normal);
  FastMat2 tangent;
};

#endif
