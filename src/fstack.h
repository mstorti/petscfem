// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: fstack.h,v 1.10 2003/02/12 19:32:09 mstorti Exp $

#ifndef FSTACK_H
#define FSTACK_H

#include <cassert>
#include <vector>
#include <string>
#include <stddef.h>
#include <regex.h>
#include <string.h>

#ifdef RH60
#include <src/libretto.h>
#endif

#include <libretto/darray.h>
#include <libretto/autostr.h>
#include <libretto/autobuf.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Allows reading from a set of files with preprocessing
    capabilities. Supports file inclusion, comments and continuation
    lines. 
    @author M. Storti */ 
class FileStack {

private: 
  /// stack of file pointers
  Darray *file_stack;
  /// the file at the top of the stack
  FILE *file_at_top,*echo_stream;
  /// Flag indicates whether lines are output to the echo stream
  int echo;
  /// the pile of unread lines
  Darray *read_buffer;
  /// The current line buffers
  Autostr *buf, *bufr,*linecopy;
  Autobuf *abufr,*abuf;
  regex_t blank_line,include;
  vector<string> file_names;
  vector<int> file_pos;
  int pos;

public:

  enum error { read_ok, cant_open, eof} last_error_m;
  error last_error();
  /** Reads a line from the file stack. 
      @author M. Storti
      @param line the line that has been read
      @return error code  */ 
  int get_line(char *& line );

  /** Gives a pointer to the last line read.
      @author M. Storti
      @param line the line that has been read  */ 
  const char * line_read(void) const;

  /** Position in file at top. 
      @return line number */ 
  int line_number() const;

  /** Name of current file at top.
      @return name of file at top */ 
  const char * file_name() const;

  /** Prints the current file stack. */ 
  void print(void) const;

  /** Unreads a line
      @author M. Storti */ 
  int unread_line(const char * line);

  /** Closes the file-stack. 
      @author M. Storti */ 
  void close();

  /** Opens a file. 
      @author M. Storti */ 
  int open(const char *filename);

  /** Constructor from the name of the main file. 
  @author M. Storti
  @param filename the name of the main file. */ 
  FileStack(const char *filename);

  /** Destructor.  
      @author M. Storti */ 
  ~FileStack(void);

  /** Check if the file has been opened correctly
      @return error code */ 
  int ok(void) {return file_at_top!=NULL;};

  /** Set stream where to write echo lines
      @param echo_s (input) the output stream */ 
  void set_echo_stream(FILE *echo_s) {echo_stream = echo_s;};
};

#endif
