//__INSERT_LICENSE__
//$Id: parsefldm.cpp,v 1.2 2001/04/14 13:20:06 mstorti Exp $

#include <cstdio>
#include <ctype.h>
#include <cmath>

#include <string>
#include <list>
#include <vector>

#include "getarray.h"
#include "parsefld.tab.h"

const char *line_to_parse;

typedef list<string> SList;
vector<SList> work_area;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void parse_fields_line(const char *line)"
void parse_fields_line(SList &fields,const char *line) {
  work_area.clear();
  line_to_parse = line;
  field_parse();
  fields.clear();
  fields.splice(fields.end(),work_area[0]);
}

void add_block(int slist_,int block_) {
  SList &slist = work_area[slist_];
  SList &block = work_area[block_];
  slist.splice(slist.end(),block);
}

void dumpwa(void) {
  for (int j=0; j<work_area.size(); j++) {
    printf("work_area[%d]: ",j);
    SList &w = work_area[j];
    for (SList::iterator k=w.begin(); k!=w.end(); k++) {
      printf("\"%s\" ",k->c_str());
    }
    printf("\n");
  }
}

int create_list(char *ident) {
  int n = work_area.size();
  work_area.resize(n+1);
  work_area[n].push_back(string(ident));
  delete[] ident;
  return n;
}

int create_ident_l(char *ident,int length) {
#define ML 20
  char nmbr[ML]; // This is ugly!!
  string s;
  int n = work_area.size();
  work_area.resize(n+1);
  for (int j=0; j<length; j++) {
    s = string(ident);
    snprintf(nmbr,ML,"[%d]",j+1);
    s.append(nmbr);
    work_area[n].push_back(s);
  }
  delete[] ident;
  return n;
}
#undef ML

int create_ident_list(char *ident,int slist_) {
  string s,b;
  SList &slist = work_area[slist_];
  b = string(ident).append(".");
  for (SList::iterator j=slist.begin(); j!=slist.end(); j++) {
    s = b;
    *j = s.append(*j);
  }
  return slist_;
}

#if 0
void add_block(void *slist_,void *block_) {
  SList &slist = *(SList *)slist_;
  SList &block = *(SList *)block_;
  for (int j=0; j<block.size(); j++) {
    slist.splice(slist.end(),block);
  }
  delete &block;
}

void *create_list(char *ident) {
  SList *slist = new SList;
  slist->push_back(string(ident));
  delete[] ident;
  return slist;
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int field_lex (void)"
int field_lex (void) {
  string s;
  int c;
  int j;
  
  c = *line_to_parse++;
  // skip whitespace 
  while (c == ' ' || c == '\t') {c = *line_to_parse++; }

  if (isalpha (c)) { 
    s="";
    j=0; 
    while (isalpha(c) || isdigit(c) || c=='_') {
      s.push_back(c);
      c = *line_to_parse++;
    }
    line_to_parse--;
    field_lval.string = new char[s.size()+1];
    strcpy(field_lval.string,s.c_str());
    // printf("creating %p with %s\n",field_lval.string,field_lval.string);
    return IDENT;
  } else if (isdigit (c)) { 
    int len;
    s="";
    j=0; 
    while (isdigit(c)) {
      s.push_back(c);
      c = *line_to_parse++;
    }
    line_to_parse--;
    sscanf(s.c_str(),"%d",&len);
    field_lval.num=len;
    return LENGTH;
  } else if (c == '[' || c == ']') {
    return c;
  } else if (c == '\n' || c == '\0') {
    return 0;
  }
  return 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ ""
int field_error (char *s) {
  printf ("parse field error: %s\n", s);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// This is a test
#if 1
#define MAXL 1000
int main() {
  char *input;
  int N=MAXL;
  input = (char *)malloc(sizeof(char)*MAXL);
  SList fields;
  
  while (1) {
    printf("enter line to parse> ");
    fgets(input,N,stdin);
    parse_fields_line(fields,input);
    int jj=0;
    for (SList::iterator j=fields.begin(); j!=fields.end(); j++) {
      printf("%d -> %s\n",jj++,j->c_str());
    }
    dumpwa();
  }
}
#endif
