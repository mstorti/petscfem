// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: syncbuff.cpp,v 1.6 2004/01/21 18:37:39 mstorti Exp $
#include <list>
#include <iostream>
#include <src/distcont.h>
#include <src/distcont2.h>
#include <algorithm>
#include <cassert>
#include <cstdio>

#include <src/syncbuff.h>
#include <src/syncbuff2.h>

using namespace std;

/// This is the stream where all elements are printed. 
FILE * KeyedLine::output = stdout;
/// Flags whether line numbers are printed. 
int KeyedLine::print_keys = 1;
/// Flags whether newlines are printed at the end of each line
int KeyedLine::print_newlines = 1;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int operator<(const KeyedLine& left, const KeyedLine& right) {
  // Simply compare the keys
  return left.key < right.key;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int KeyedLine::size_of_pack() const {
  // To pack, we could store only the key and the
  // string, since we could deduce the length in unpack
  // with strlen. But we add the length for better checking. 
  return 2*sizeof(int)+strlen(line)+1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedLine::pack(char *&buff) const {
  // Pack the key
  memcpy(buff,&key,sizeof(int));
  buff += sizeof(int);

  // Pack the string length
  int len = strlen(line);
  memcpy(buff,&len,sizeof(int));
  buff += sizeof(int);

  // Pack the string itself
  memcpy(buff,line,strlen(line)+1);
  buff += strlen(line)+1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedLine::unpack(const char *& buff) {
  // Unpack the key
  memcpy(&key,buff,sizeof(int));
  buff += sizeof(int);

  // Unpack the string length
  int len;
  memcpy(&len,buff,sizeof(int));
  buff += sizeof(int);

  // Unpack the string
  assert(!line);
  line = new char[len+1];
  memcpy(line,buff,len+1);
  buff += strlen(line)+1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedLine::print() {
  // Prints the line with or without the keys. 
  if (print_keys) 
    fprintf(output,"%d: %s",key,line);
  else 
    fprintf(output,"%s",line);
  if (print_newlines) fprintf(output,"\n");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
KeyedLine::KeyedLine(const KeyedLine &kl) {
  // Copys constructor
  // Cosnider the self copy case. 
  if (this==&kl) return;
  key = kl.key;
  // Consider the special case of a void object
  if (kl.line) {
    line = new char[strlen(kl.line)+1];
    memcpy(line,kl.line,strlen(kl.line)+1);
  } else line=NULL;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedLine::build(int k,const char *s) {
  // Internal function that builds an instance from key and C-string. 
  key = k;
  int len = strlen(s);
  line = new char[len+1];
  memcpy(line,s,len+1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Ctor 
KeyedLine::KeyedLine(int k,const char *s) { build(k,s); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Ctor 
KeyedLine::KeyedLine(int k,const AutoString &as) { 
  build(k,as.str()); 
}

// This is the tricky part due to the problem
// with partial template specializations
SYNC_BUFFER_FUNCTIONS(KeyedLine);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Special shortcut to the push functions
void KeyedOutputBuffer::push(int k,const AutoString &as) {
  push_back(KeyedLine(k,as.str()));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Special shortcut to the push functions
void KeyedOutputBuffer::push(int k,const char *s) {
  push_back(KeyedLine(k,s));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedOutputBuffer::flush() {
  // Print and clear
  print();
  clear();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedOutputBuffer::vprintf(const char * tmplt,va_list ap) { 
  as.vsprintf(tmplt,ap); 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedOutputBuffer::printf(const char * tmplt, ...) { 
  va_list ap;
  va_start(ap,tmplt);
  vprintf(tmplt,ap);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedOutputBuffer::vcat_printf(const char * tmplt,va_list ap) { 
  as.vcat_sprintf(tmplt,ap); 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedOutputBuffer::cat_printf(const char * tmplt, ...) { 
  va_list ap;
  va_start(ap,tmplt);
  vcat_printf(tmplt,ap);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void KeyedOutputBuffer::push(int k) { 
  push(k,as); as.clear(); 
}
