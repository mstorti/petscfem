// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: syncbuff.cpp,v 1.2 2004/01/11 16:22:04 mstorti Exp $
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
void KeyedLine::build(int k,const char *s) {
  key = k;
  int len = strlen(s);
  line = new char[len+1];
  memcpy(line,s,len+1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
KeyedLine::KeyedLine(int k,const char *s) { build(k,s); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
KeyedLine::KeyedLine(int k,const AutoString &as) { 
  build(k,as.str()); 
}

SYNC_BUFFER_FUNCTIONS(KeyedLine);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedOutputBuffer::push(int k,const AutoString &as) {
  push_back(KeyedLine(k,as.str()));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedOutputBuffer::push(int k,const char *s) {
  push_back(KeyedLine(k,s));
}
