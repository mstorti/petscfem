// -*- mode: C++ -*-

/*
  This file belongs to the PETSc - FEM package, a library and
  application suite oriented to the Finite Element Method based on PETSc. 
  Copyright (C) 1999, 2000  Mario Alberto Storti
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.

*/

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
