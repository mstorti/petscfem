// -*- mode: c++ -*-
/*__INSERT_LICENSE__*/
//$Id: util2.h,v 1.8 2001/06/23 16:43:14 mstorti Exp $
#ifndef UTIL2_H
#define UTIL2_H

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
#if 0
  /// create a chronometer
  Chrono();
  /// destroy chronometer
  ~Chrono() {};
#endif
  /// return elapsed CPU time from start
  double elapsed() const {return ((double) (clock() - start_time)) / CLOCKS_PER_SEC;};
  /// reset start time to actual time
  void start() {start_time = clock();};
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

#endif
