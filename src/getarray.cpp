//__INSERT_LICENSE__
//$Id: getarray.cpp,v 1.9 2003/07/02 03:36:13 mstorti Exp $

#include <cstdio>
#include <ctype.h>
#include <cmath>
#include <cstring>

#include <string>
#include <vector>

#include "getarray.h"
#include "getarrgr.tab.h"

using namespace std;

vector<int> *prop_len;
vector<string> *prop_name;
const char *line_to_parse;

void add_entry(char *s,int l) {
  prop_len->push_back(l);
  prop_name->push_back(string(s));
  // printf("deleting %p with %s in add_entry\n",s,s);
  delete[] s;
}

int yylex (void) {
  string s;
  int c;
  
  c = *line_to_parse++;
  // skip whitespace 
  while (c == ' ' || c == '\t') {c = *line_to_parse++; }

  if (isalpha (c)) { 
    s="";
    while (isalpha(c) || isdigit(c) || c=='_') {
      s += c;
      c = *line_to_parse++;
    }
    line_to_parse--;
    yylval.string = new char[s.size()+1];
    strcpy(yylval.string,s.c_str());
    // printf("creating %p with %s\n",yylval.string,yylval.string);
    return IDENT;
  } else if (isdigit (c)) { 
    int len;
    s="";
    while (isdigit(c)) {
      s += c;
      c = *line_to_parse++;
    }
    line_to_parse--;
    sscanf(s.c_str(),"%d",&len);
    yylval.num=len;
    return LENGTH;
  } else if (c == '[' || c == ']') {
    return c;
  } else if (c == '\n' || c == '\0') {
    return 0;
  }
  return 1;
}

void parse_props_line(const char *line,vector<string> &prop_name_,vector<int>
		      &prop_len_) {
  prop_name = &prop_name_;
  prop_len = &prop_len_;
  prop_name->clear();
  prop_len->clear();
  line_to_parse = line;
  yyparse();
}

/* Called by yyparse on error */
int yyerror (char *s) {
  printf ("%s\n", s);
  return 0;
}

//void opt_def_sl_action

#if 0
void add_index(void **sl,char *s,int first) {
  vector<string> *sll;
  if (first) {
    sll = new vector<string>;
    *sl = (void *) sll;
  } else {
    sll = (vector<string> *)*sl;
  }
  //printf("adding to pointer %p, string %s, first %d\n",sll,s,first);
  sll->push_back(string(s));
  // printf("deleting %p with %s\n",s,s);
  delete[] s;
}

void opt_def_sl_action(char *s,void *sl) {
  // printf("option with string indices -> %s[%p]\n",s,sl);
  assert(0); // not string subscripts yet
  vector<string> *sll = (vector<string> *) sl;
//    for (int j=0; j<sll->size(); j++) {
//      printf("%s.%s\n",s,(*sll)[j].c_str());
//    }
  delete sll;
  delete[] s;
}
#endif

// This is a test
#if 0

double drand() {  
  return ((double)(rand()))/((double)(RAND_MAX));
}

int irand(int imin,int imax) {
  return int(rint(drand()*double(imax-imin+1)-0.5))+imin;
}

int main() {
  char input[1000];
  vector<string> names;
  vector<int> lens;
  while (1) {
    sprintf(input,"   pp%d   [  qq%d   rr%d  eeee%d   fffff%d   ]   ",
	    irand(1,100),irand(1,100),irand(1,100),irand(1,100),irand(1,100));
    printf("%s\n",input);
    parse_props_line(input,names,lens);
    for (int j=0; j<names.size(); j++) {
      printf("name %s, len %d\n",names[j].c_str(),lens[j]);
    }
  }
}
#endif
