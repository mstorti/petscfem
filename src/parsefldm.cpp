//__INSERT_LICENSE__
//$Id: parsefldm.cpp,v 1.5 2001/04/15 22:39:03 mstorti Exp $

#include <cstdio>
#include <ctype.h>
#include <cmath>

#include <string>
#include <list>
#include <vector>

#include "wspace.h"
#include "getarray.h"
#include "nodedata.h"
#include "parsefld.tab.h"

const char *line_to_parse;

WorkSpace<SList> work_area;
WorkSpace<string> string_space;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void dumpwa(int debug=0)"
void dumpwa(int debug=0) {
  if (debug) printf("work_area: ----------\n");
  int jj=0;
  for (WorkSpace<SList>::iterator j=work_area.begin(); j!=work_area.end(); j++) {
    SList &w = *(SList *)(*j);
    if (w.size()==0 && !debug) continue;
    printf("work_area[%d]: ",++jj);
    for (SList::iterator k=w.begin(); k!=w.end(); k++) {
      printf("\"%s\" ",k->c_str());
    }
    printf("\n");
  }
  if (debug) {
    printf("string_space: ----------\n");
    int jj=0;
    for (WorkSpace<string>::iterator j=string_space.begin(); 
	 j!=string_space.end(); j++) {
      printf("%d -> \"%s\"\n",++jj,(*j)->c_str());
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void parse_fields_line(const char *line)"
int parse_fields_line(SList &fields,const char *line) {
  work_area.clear();
  line_to_parse = line;
  int ierr = field_parse();
  if (ierr>0) return ierr;
  fields.clear();
  if (work_area.size()>1) {
    printf("error: work_area.size(): %d\n",work_area.size());
    abort();
  }
  if (string_space.size()) {
    printf("error: string_space.size(): %d\n",);
    abort();
  }
  fields.splice(fields.end(),**work_area.begin());
  dumpwa(0);
  work_area.clear();
  string_space.clear();

  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void add_block(void *slist_,void *block_)"
void *add_block(void *slist_,void *block_) {
  SList &slist = *(SList *)slist_;
  SList &block = *(SList *)block_;
  slist.splice(slist.end(),block);
  work_area.free(&block);
  return &slist;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int create_list(void *ident_)"
void *create_list(void *ident_) {
  string &ident = *(string *)ident_;
  SList &sl = *work_area.newobj();
  sl.push_back(ident);
  string_space.free(&ident);
  return &sl;
}

void *f_id_subs_l(void *ident_,void *slist_) {
  string &ident = *(string *)ident_;
  SList &slist = *(SList *)slist_;
  string s;
  for (SList::iterator j=slist.begin(); j!=slist.end(); j++) {
    s=ident;
    s.append(*j);
    *j = s;
  }
  string_space.free(&ident);
  return slist_;
}
  
void *list_product(void *sl1_,void *sl2_) {
  SList &sl1 = *(SList *)sl1_;
  SList &sl2 = *(SList *)sl2_;
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
  work_area.free(&sl2);
  return &sl1;
}
      
void *dotify(void *sl_) {
  SList &sl = *(SList *)sl_;
  string s;
  for (SList::iterator j=sl.begin(); j!=sl.end(); j++) {
    s= string(".");
    s.append(*j);
    *j = s;
  }
  return &sl;
}
    
void *create_list_from_length(int length) {
  SList &sl = *work_area.newobj();
  
  string s;
#define ML 20
  char nmbr[ML]; // This is ugly!!

  for (int j=0; j<length; j++) {
    snprintf(nmbr,ML,"[%d]",j+1);
    s = string(nmbr);
    sl.push_back(s);
  }
  return &sl;
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
    string *ss = string_space.newobj();
    *ss = s;
    field_lval.gen_ptr = ss;
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

#else

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// This is a test
int main() {
  WorkSpace<double> wd;
  double *a = wd.newobj();
  double *b = wd.newobj();
  double *c = wd.newobj();
  
  wd.free(a);
  wd.free(b);
  wd.free(c);
  
}
#endif
