// -*- mode: C++ -*-
//__INSERT_LICENSE__
// $Id: autostr.h,v 1.7 2003/05/04 16:34:54 mstorti Exp $
#ifndef PETSCFEM_AUTOSTR_H
#define PETSCFEM_AUTOSTR_H
#include <cstdarg>
#include <cstdio>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** A string class that dynamically grows, and allows safe use of the
      asprintf function, and friends. The #asprintf# is used instead of
      the standard #sprintf# function so that overflow of the array is
      avoided. In addition all functions are reentrant. Many functions 
      are variadic (see doc on variadic functions in the GNU glibc 
      documentation). For many functions there is the #foo(...)# and
      #vfoo(va_list ap)# version of the function. The variable #tmplt# 
      is a standard #printf# template for output conversion (with #%s# 
      etc...)
*/
class AutoString {
private:
  /// The actual string
  char *s;
  /// The number of bytes allocated in #s#
  int n;
  /// Total number of characters allocated by the class
  static int total;
public:
  /// Ctor.
  AutoString();
  /// Dtor. The string is initialized to the empty string #""#. 
  ~AutoString();
  /** Allows const access to the C buffer. 
      @return const pointer to internal buffer*/ 
  const char *str() const;
  /** Returns size of buffer (may be not equal to the length of 
      the string actually stored. 
      @return size of internal buffer */
  int size() const;
  /** Size of string as detected by the position of the #\0# terminator. 
      @return length of string */ 
  int len() const;
  /** Resizes the internal buffer. 
      @param m (input) new size */ 
  void resize(int m);
  /** Sets to the null string #""#. 
      @return reference to self */ 
  AutoString & clear();
  /** Sets string to result of #sprintf#. #as.sprintf(patt,...)# 
      is equivalent to #sprintf(as,patt,...)# if #as# were a normal
      C string. 
      @param tmplt (input) output template
      @return pointer to self */ 
  AutoString & sprintf(const char * tmplt,...);
  /** Variadic explicit version of #sprintf(const char * tmplt,...)#. 
      @param NAME (input) TEXT
      @return pointer to self */ 
  AutoString & vsprintf(const char * tmplt,va_list ap);
  /** Appends to string the result of the #sprintf# call. 
      @param tmplt (input) output template 
      @return pointer to self */ 
  AutoString & cat_sprintf(const char *,...);
  /** Appends to string the result of the #sprintf# call. 
      Variadic explicit version of #cat_sprintf()#. 
      @param tmplt (input) output template 
      @param ap (input) array of arguments 
      @return pointer to self */ 
  AutoString & vcat_sprintf(const char *,va_list ap);
  /** Appends #s# to string. 
      @param s (input) string to be appended
      @return pointer to self */ 
  AutoString & cat(const AutoString &s);
  /** Appends #s# to string. 
      @param s (input) string to be appended
      @return pointer to self */ 
  AutoString & cat(const char *s);
  /** Prints to file. 
      @param fid (input) file identifier
      @return pointer to self */ 
  AutoString & fprintf(FILE *fid,...);
  /** Prints to file. Variadic explicit version of #fprintf()#. 
      @param fid (input) file identifier
      @return pointer to self */ 
  AutoString & vfprintf(FILE *fid,va_list ap);
  /** Prints to standard output. 
      @return pointer to self */ 
  AutoString & print();
  /** Prints to standard output. (const version). 
      @return pointer to self */ 
  const AutoString & print() const;
  /** Sets string to #s#. 
      @param s (input) string to be copied to. 
      @return pointer to self */ 
  AutoString & set(const AutoString &s);
  /** Sets string to #s#. 
      @param s (input) string to be copied to. 
      @return pointer to self */ 
  AutoString & set(const char *s);
  /** Reads a line from file #fid#. 
      @param fid (input) file for reading. 
      @return number of characters read. */ 
  int getline(FILE *fid);
};

#endif
