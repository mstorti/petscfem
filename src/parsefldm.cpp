//__INSERT_LICENSE__
//$Id: parsefldm.cpp,v 1.3 2001/04/14 21:00:48 mstorti Exp $

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
#define __FUNC__ "void dumpwa(int debug=0)"
void dumpwa(int debug=0) {
  for (int j=0; j<work_area.size(); j++) {
    SList &w = work_area[j];
    if (w.size()==0 && !debug) continue;
    printf("work_area[%d]: ",j);
    for (SList::iterator k=w.begin(); k!=w.end(); k++) {
      printf("\"%s\" ",k->c_str());
    }
    printf("\n");
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void parse_fields_line(const char *line)"
void parse_fields_line(SList &fields,const char *line) {
  work_area.clear();
  line_to_parse = line;
  field_parse();
  fields.clear();
  fields.splice(fields.end(),work_area[0]);
  dumpwa(0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void add_block(int slist_,int block_)"
void add_block(int slist_,int block_) {
  SList &slist = work_area[slist_];
  SList &block = work_area[block_];
  slist.splice(slist.end(),block);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int create_list(char *ident)"
int create_list(char *ident) {
  int n = work_area.size();
  work_area.resize(n+1);
  work_area[n].push_back(string(ident));
  delete[] ident;
  return n;
}

int f_id_subs_l(char *ident,int slist_) {
  SList &slist = work_area[slist_];
  string s,b = string(ident);
  for (SList::iterator j=slist.begin(); j!=slist.end(); j++) {
    s=b;
    s.append(*j);
    *j = s;
  }
  return slist_;
}
  
int list_product(int sl1_,int sl2_) {
  SList &sl1 = work_area[sl1_];
  SList &sl2 = work_area[sl2_];
  SList result;
  result.clear();
  string s;
  for (SList::iterator j1=sl1.begin(); j1!=sl1.end(); j1++) {
    for (SList::iterator j2=sl2.begin(); j2!=sl2.end(); j2++) {
      s=(*j1);
      s.append(*j2);
      result.push_back(s);
    }
  }
  sl1.clear();
  sl1.splice(sl1.end(),result);
  sl2.clear();
  return sl1_;
}
      
int dotify(int sl_) {
  SList &sl = work_area[sl_];
  string s;
  for (SList::iterator j=sl.begin(); j!=sl.end(); j++) {
    s= string(".");
    s.append(*j);
    *j = s;
  }
  return sl_;
}
    
int create_list_from_length(int length) {
  int n = work_area.size();
  work_area.resize(n+1);
  SList &sl = work_area[n];
  
  string s;
#define ML 20
  char nmbr[ML]; // This is ugly!!

  for (int j=0; j<length; j++) {
    snprintf(nmbr,ML,"[%d]",j+1);
    s = string(nmbr);
    sl.push_back(s);
  }
  return n;
}
#undef ML

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
#undef __FUNC__
#define __FUNC__ "int create_ident_l(char *ident,int length)"
int make_length_list(int length,int list_=-1) {
#define ML 20
  char nmbr[ML]; // This is ugly!!
  string s;
  int n = work_area.size();
  work_area.resize(n+1);
  for (int j=0; j<length; j++) {
    snprintf(nmbr,ML,"[%d]",j+1);
    s = string(nmbr);
    work_area[n].push_back(s);
  }
  return list_;
}
#undef ML
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int create_ident_list(char *ident,int slist_)"
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
  }
}
#endif
