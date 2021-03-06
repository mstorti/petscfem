// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: syncbuff.h,v 1.12 2007/01/30 19:03:44 mstorti Exp $
#include <list>
#include <iostream>
#include <src/distcont.h>
#include <src/distcont2.h>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <src/autostr.h>

extern int SIZE, MY_RANK;

using namespace std;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Trivial partitioning class for use with the 
    syncronized buffer classes. All elements belong 
    to processor 0. 
    @param NAME (input) TEXT
    @return RETURN-TEXT */ 
template<typename T>
class TrivialPartitioner {
public:
  /** All elements belong to processor 0. 
      @param t (input) the element to be scattered
      @param nproc (output) number of processors to 
      which this element belongs
      @param plist (output) user sets ths first nproc values 
      of this array to the number of processes owners */
  void processor(const T &t,int &nproc,int *plist);
};

/** I had trouble in make a partial specialization of #DistCont<>#
    for #DistCont<list<T>,T,TrivialPartitioner<T> ># 
    and I finally wrote this macro for defining the wrappers to 
    the pack/unpack functions specific to this case. 
    @param T (input) the ValueType class
    @return definitions of functions */ 
#define SB(T) DistCont<list<T>,T,TrivialPartitioner<T> >

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** A distributed container based on #list# that receives elements
    and then a #scatter# sends all of them to the master node, 
    probably for printing.  */ 
template<typename T>
class SyncBuffer : public DistCont<list<T>,T,TrivialPartitioner<T> > {
  /// The partitioner (all elements belong to the master)
  TrivialPartitioner<T> part;
public:
  /** Flags whether the lines should be sorted by key in the 
      output */
  int sort_by_key;
  /// Constructor
  SyncBuffer() : 
    DistCont<list<T>,T,TrivialPartitioner<T> >(&part),
    sort_by_key(1) {}
  /// Prints all elements to the output stream
  void print();
  /// This is for checking the pack/unpack routines. 
  void check_pack();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define SYNC_BUFFER_FUNCTIONS(T)		\
						\
template<> int SB(T)				\
::size_of_pack(const T &t) const {		\
  return t.size_of_pack();			\
}						\
						\
template<> void SB(T)				\
::pack(const T &t, char *&buff) const {		\
  t.pack(buff);					\
}						\
						\
template<> void SB(T)				\
::unpack(T &t,const char *& buff) {		\
  t.unpack(buff); }				\
						\
template<> void SB(T)				\
::combine(const T &t) {				\
push_back(t); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** A line with an integer key (see doc for 
    #KeyedOutputBuffer#) and a string. All lines are ordered by the
    key and printed to the output stream. */ 
class KeyedLine {
private:
  /// The key to order lines
  int key;
  /// The string to be printed
  char *line;
  /** Builds an instance from an integer key and a external string
      @param key (input) the key of the line
      @param s (input) the line of the text */ 
  void build(int key,const char *s);
public:
  /// Flags whether the keys are printed to the output
  static int print_keys;
  /// Flags whether newlines are printed at the end of each line
  static int print_newlines;
  /// The output stream
  static FILE *output;
  /// Default constructor 
  KeyedLine() : key(0), line(NULL) {}
  /// Copy Ctor
  KeyedLine(const KeyedLine &kl);
  /// Constructor from autostring
  KeyedLine(int key,const AutoString &as);
  /// Constructor from key and C-string. 
  KeyedLine(int key,const char *s);
  /// Dtor.
  ~KeyedLine() { if (line) delete[] line; }
  /// Used for sorting the lines
  friend int operator<(const KeyedLine& left, const KeyedLine& right);
  /**  Call back for the distributed container. Returns
       size of the buffer needed to store this element.  */ 
  int size_of_pack() const;
  /** Effectively packs the object into the buffer, 
      upgrading the pointer. 
      @param buff (input/output) the output buffer. 
      @return RETURN-TEXT */ 
  void pack(char *&buff) const;
  /** Extracts the object from the buffer, upgading the buffer. 
      @param buff (input) the output buffer */ 
  void unpack(const char *& buff);
  /// Print the object
  void print();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Synchronized printing buffer. Objects are added with 
    a #push# operation and then a #print# operation does the 
    scatter, sorting according to the keys and printing.  */ 
class KeyedOutputBuffer : public SyncBuffer<KeyedLine> {
private:
  AutoString as;
public:
  /** Sets buffer to string created like #printf#. 
      @param tmplt (input) the printf-like template 
      @param VARARGS (input) the items to be printed. */
  void printf(const char * tmplt, ...);
  /** Explicit variadic form of #printf(...)#
      @param tmplt (input) the printf-like template 
      @param ap (input) variadic vector of items to be printed. */
  void vprintf(const char * tmplt,va_list ap);
  /** Appends to buffer a string created like #printf#. 
      @param tmplt (input) the printf-like template 
      @param VARARGS (input) the items to be printed. */
  void cat_printf(const char * tmplt, ...);
  /** Explicit variadic form of #vcat_printf(...)#
      @param tmplt (input) the printf-like template 
      @param ap (input) variadic vector of items to be printed. */
  void vcat_printf(const char * tmplt,va_list ap);
  /** Flushes the current string in the internal buffer to the
      buffer with given key. 
      @param key (input) the key associated with the current string. */
  void push(int k);
  /** Inserts a line in the buffer. 
      @param k (input) the key of the line
      @param s (input) the line of text */ 
  void push(int k,const AutoString &s);
  /** Inserts a line in the buffer. 
      @param k (input) the key of the line
      @param s (input) the line of text */ 
  void push(int k,const char *s);
  /// Makes the scatter and prints the output. 
  void flush();
};
