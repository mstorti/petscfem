// -*- mode: C++ -*-
//__INSERT_LICENSE__
//$Id: fstack.h,v 1.3 2001/04/01 01:35:06 mstorti Exp $

#ifndef FSTACK_H
#define FSTACK_H

#include <cassert>
#include <stddef.h>
#include <regex.h>
#include <string.h>

#ifdef RH60
#include "libretto.h"
#endif

#include <libretto/darray.h>
#include <libretto/autostr.h>
#include <libretto/autobuf.h>

//#include <petsc.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Allows reading from a set of files with preprocessing
    capabilities. Supports file inclusion, comments and continuation
    lines. 

    @author M. Storti
*/ 
class FileStack {

private: 
  /// stack of file pointers
  Darray *file_stack;
  /// the file at the top of the stack
  FILE *file_at_top;
  /// the pile of unread lines
  Darray *read_buffer;
  /// The current line buffers
  Autostr *buf, *bufr;
  Autobuf *abufr,*abuf;
  regex_t blank_line,include;

public:
  /** Reads a line from the file stack. 
      @author M. Storti
      @param line the line that has been read
  */ 
  int get_line(char *& line );

  /** Unreads a line
      @author M. Storti
  */ 
  int unread_line(const char * line);

  /** Closes the file-stack. 
      @author M. Storti
  */ 
  void close();

  /** Opens a file. 
      @author M. Storti
  */ 
  void open(const char *filename);

  /** Constructor from the name of the main file. 
  @author M. Storti
  @param filename the name of the main file. 
  */ 
  FileStack(const char *filename);

  /** Destructor.  
      @author M. Storti
  */ 
  ~FileStack(void);

};

#endif
