// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: syncbuff.cpp,v 1.1 2004/01/11 15:58:46 mstorti Exp $
#include <list>
#include <iostream>
#include <src/distcont.h>
#include <src/distcont2.h>
#include <src/debug.h>
#include <algorithm>
#include <cassert>
#include <cstdio>

#include <src/syncbuff.h>

using namespace std;

FILE * KeyedLine::output = stdout;
int KeyedLine::print_line_numbers = 1;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int operator<(const KeyedLine& left, const KeyedLine& right) {
  return left.key < right.key;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int KeyedLine::size_of_pack() const {
  return 2*sizeof(int)+strlen(line)+1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedLine::pack(char *&buff) const {
  memcpy(buff,&key,sizeof(int));
  buff += sizeof(int);

  int len = strlen(line);
  memcpy(buff,&len,sizeof(int));
  buff += sizeof(int);

  memcpy(buff,line,strlen(line)+1);
  buff += strlen(line)+1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedLine::unpack(const char *& buff) {
  memcpy(&key,buff,sizeof(int));
  buff += sizeof(int);

  int len;
  memcpy(&len,buff,sizeof(int));
  buff += sizeof(int);
  assert(!line);
  line = new char[len+1];
  memcpy(line,buff,len+1);
  buff += strlen(line)+1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedLine::print() {
  if (print_line_numbers) 
    fprintf(output,"%d: %s\n",key,line);
  else 
    fprintf(output,"%s\n",line);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
KeyedLine::KeyedLine(const KeyedLine &kl) {
  if (this==&kl) return;
  key = kl.key;
  if (kl.line) {
    line = new char[strlen(kl.line)+1];
    memcpy(line,kl.line,strlen(kl.line)+1);
  } else line=NULL;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
KeyedLine::KeyedLine(int k,const AutoString &as) {
  key = k;
  const char *l = as.str();
  int len = strlen(l);
  line = new char[len+1];
  memcpy(line,l,len+1);
}

SYNC_BUFFER_FUNCTIONS(KeyedLine);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedOutputBuffer::push(int k,const AutoString &s) {
  push_back(KeyedLine(k,s));
}

