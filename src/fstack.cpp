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
  abuf_destroy(abuf);
  abuf_destroy(abufr);
  da_destroy(read_buffer);
  da_destroy(file_stack);
  regfree(&blank_line);
  regfree(&include);
}

// Temporally. I write close() equal to the destructor, but in the
// future it will close all the files in the stack. 
#undef __FUNC__
#define __FUNC__ "FileStack::close" 
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FileStack::close(void) {
  if (da_length(file_stack)>0 || da_length(read_buffer)>0) {
    printf("FileStack::close(): FileStack not empty!!\n"
	   "%d files open\n"
	   "%d lines in the buffer\n",
	   da_length(file_stack),da_length(read_buffer));
  }
  fclose(file_at_top);
  while (da_length(file_stack)>0) {
    da_pop(file_stack,&file_at_top);
    fclose(file_at_top);
  }
}

#undef __FUNC__
#define __FUNC__ "FileStack::FileStack" 
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FileStack::FileStack(const char *filename) {

  // char line[LINESIZE];
  file_at_top = fopen(filename,"r");
  if (!file_at_top) {
    printf("Couldn't open file \"%s\"!!\n",filename);
    exit(0);
  }
  file_stack = da_create(sizeof(FILE *));
  read_buffer = da_create(sizeof(Autobuf *));
  buf = astr_create();
  bufr = astr_create();
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
void FileStack::open(const char *filename) {
  // char line[LINESIZE];
  file_at_top = fopen(filename,"r");
  if (!file_at_top) {
    printf("Couldn't open file \"%s\"!!\n",filename);
    exit(0);
  }
}

#undef __FUNC__
#define __FUNC__ "FileStack::get_line" 
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int FileStack::get_line(const char *& line) {
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
    return 0;
  }

  int nread;
  while(1) {
    nread = astr_getline(bufr,file_at_top);

    if (nread<0) {
      // If couldn't read then pop the following file in the stack
      if (da_length(file_stack)>0) {
	fclose(file_at_top);
	da_pop(file_stack,&file_at_top);
	continue;
      } else {
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
      file_at_top = fopen(token,"r");

      if (!file_at_top) {
	printf("Couldn't open file \"%s\"!!\n",token);
	exit(0);
      }
      
      continue;
    }
      
    readlen = strlen(bufrp);
    if (readlen>1) {
      contnd=0;
      if (bufrp[readlen-1] == '\\' ) {
	readlen -= 1;
	contnd=1;
	bufrp[readlen-1] = '\0';
      }
      lenn = len+ readlen-1;
//        if (lenn > MAX_LINE) {
//  	printf("exceeded size of line: %d\n",MAX_LINE);
//  	exit(1);
//        }

      abuf_cat_s (abuf,bufrp);
      len=lenn;
    }
    if ( ! contnd ) break;
  }
  //  abuf_data(abuf)[len+1]='\0';
  abuf_cat_c (abuf,'\0');
  line = (char *)(abuf_data(abuf));
  ///  printf("read: ->\"%s\"\n",line);
  return 0;
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
}
  
