//__INSERT_LICENSE__
//$Id: fstack.cpp,v 1.16 2004/09/25 23:11:39 mstorti Exp $
#include <stdlib.h>
#include "fstack.h"

// for regexp delimiter
#define BSP "[ \t\n]"

#undef __FUNC__
#define __FUNC__ "FileStack::~FileStack" 
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FileStack::~FileStack(void) {
  close();
  astr_destroy(buf);
  astr_destroy(bufr);
  astr_destroy(linecopy);
  abuf_destroy(abuf);
  abuf_destroy(abufr);
  da_destroy(read_buffer);
  da_destroy(file_stack);
  regfree(&blank_line);
  regfree(&include);
  file_names.clear();
  file_pos.clear();
}

// Temporally. I write close() equal to the destructor, but in the
// future it will close all the files in the stack. 
#undef __FUNC__
#define __FUNC__ "FileStack::close" 
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FileStack::close(void) {
  if (da_length(file_stack)>0 || da_length(read_buffer)>0) {
    if (!quiet)
      printf("FileStack::close(): FileStack not empty!!\n"
	     "%d files open\n"
	     "%d lines in the buffer\n",
	     da_length(file_stack),da_length(read_buffer));
  }
  if (file_at_top) fclose(file_at_top);
  file_at_top = NULL;
  while (da_length(file_stack)>0) {
    da_pop(file_stack,&file_at_top);
    fclose(file_at_top);
  }
}

#undef __FUNC__
#define __FUNC__ "FileStack::FileStack" 
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FileStack::FileStack(const char *filename) : 
  quiet(0),
  echo_stream (NULL), 
  echo(0), 
  file_at_top(NULL),
  last_error_m(read_ok) {

  open(filename);
  if (!ok()) {
    last_error_m = cant_open;
    return;
  }
  file_stack = da_create(sizeof(FILE *));
  read_buffer = da_create(sizeof(Autobuf *));
  buf = astr_create();
  bufr = astr_create();
  linecopy = astr_create();
  abuf = abuf_create();
  abufr = abuf_create();
  
  int ierr = regcomp(&blank_line,"^[ \t\n]*$",REG_NOSUB);
  assert(ierr==0);
  ierr = regcomp(&include,"^__INCLUDE__",REG_NOSUB);
  assert(ierr==0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void FileStack::open(const char *)"
int FileStack::open(const char *filename) {
  // char line[LINESIZE];
  file_at_top = fopen(filename,"r");
  if (!file_at_top) {
    if (!quiet) 
      printf("Couldn't open file \"%s\"!!\n",filename);
    last_error_m = cant_open;
    return 1;
  }
  file_names.push_back(string(filename));
  file_pos.push_back(pos);
  pos=0;
  last_error_m = read_ok;
  return 0;
}

#undef __FUNC__
#define __FUNC__ "FileStack::get_line" 
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int FileStack::get_line(char * & line) {
  char ch, *chp, *token, *bufrp;
  int len,lenn,readlen,contnd,ierr;
  len=0;

  // resets line
  //if (abuf_length(abuf)>0) abuf_delete (abuf,0,-1);
  abuf_zero(abuf);

  if (da_length(read_buffer)>0) {
    abuf_destroy(abuf);
    da_pop(read_buffer,&abuf);
    line = (char *)(abuf_data(abuf));
    last_error_m = read_ok;
    return 0;
  }

  int nread;
  while(1) {
    nread = astr_getline(bufr,file_at_top);
    // This is for DOS files :-(
    if (nread>0 && astr_chars(bufr)[nread-1]=='\r') {
      astr_delete(bufr,nread-1,1);
      nread--;
    }
    pos++;

    if (nread<0) {
      // If couldn't read then pop the following file in the stack
      if (da_length(file_stack)>0) {
	fclose(file_at_top);
	da_pop(file_stack,&file_at_top);
	file_names.pop_back();
	pos = file_pos.back();
	file_pos.pop_back();
	continue;
      } else {
	last_error_m = eof;
	return 1;
      }
    }

    if (astr_length(bufr)==0) continue;
    int indx =  astr_find_c(bufr,0,'#');
    if (indx!= -1) {
      astr_delete (bufr,indx,-1);
    }

    //if (abuf_length(abufr)>0) abuf_delete (abufr,0,-1);
    //abuf_zero(abufr);
    abuf_copy_s (abufr,astr_chars(bufr));
    abuf_cat_c (abufr,'\0');
    bufrp = (char *)abuf_data(abufr);

    // flag echo
    if (!strcmp("__ECHO_ON__",bufrp)) {
      echo=1;
      continue;
    }

    // flag no echo
    if (!strcmp("__ECHO_OFF__",bufrp)) {
      echo=0;
      continue;
    }

    if (echo && echo_stream) {
      fprintf(echo_stream,"%s\n",bufrp);
    }

    // skip comments
    if (bufrp[0] == '#') continue;

    // skip blank lines
    if (!regexec (&blank_line,bufrp,0,NULL,0)) continue;

    // include file
    if (!regexec (&include,bufrp,0,NULL,0)) {
      da_push(file_stack,&file_at_top);

      // read filename
      token = strtok((char *)abuf_data(abufr),BSP);
      token = strtok(NULL,BSP); 
      string token_cpy(token);
      if (token_cpy[0]=='"') {
	int m = strlen(token);
	token_cpy.erase(m-1,1);
	token_cpy.erase(0,1);
      }
      file_at_top = fopen(token_cpy.c_str(),"r");
	
      if (!file_at_top) {
	if (!quiet)
	  printf("Couldn't open file \"%s\"!!\n",token);
	last_error_m = cant_open;
	return 1;
      }
      file_names.push_back(token_cpy);
      file_pos.push_back(pos);
      pos=0;
      
      continue;
    }
      
    readlen = strlen(bufrp);
    if (readlen>0) {
      contnd=0;
      if (bufrp[readlen-1] == '\\' ) {
	readlen -= 1;
	contnd=1;
	bufrp[readlen-1] = '\0';
      }
      lenn = len+ readlen-1;

      abuf_cat_s (abuf,bufrp);
      len=lenn;
    }
    if ( ! contnd ) break;
  }
  abuf_cat_c (abuf,'\0');
  line = (char *)(abuf_data(abuf));
  astr_copy_s(linecopy, line);
  
  last_error_m = read_ok;

  return 0;
} 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int FileStack::line_read(const char *& line )"
const char * FileStack::line_read(void) const {
  return astr_chars(linecopy);
  return 0;
}

void FileStack::print(void) const {
  if (!quiet) {
    printf("File stack:\n");
    int nfiles= file_names.size();
    for (int j=0; j<file_names.size(); j++) {
      printf("pos %d in stack: %s:%d\n",
	     j,file_names[j].c_str(),( j != nfiles-1 ? file_pos[j] : pos));
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ ""
int FileStack::unread_line(const char * line) {
  Autobuf *linebuf;
  linebuf = abuf_create();
  abuf_copy_s (linebuf,line);
  abuf_cat_c (linebuf,'\0');
  da_push(read_buffer,&linebuf);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:   
int FileStack::line_number() const {return pos;};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
const char * FileStack::file_name() const {
  return file_names.back().c_str();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FileStack::error FileStack::last_error() { return last_error_m; }
