//__INSERT_LICENSE__
//$Id: tryme.cpp,v 1.3 2003/01/08 15:54:25 mstorti Exp $
#include <stdio.h>
#include <ctype.h>

#include <string>
#include <vector>

#include "getarray.h"
#include "getarraypars.tab.h"

vector<int> prop_len;
vector<string> prop_name;

void add_entry(char *s,int l) {
  prop_len.push_back(l);
  prop_name.push_back(string(s));
}

int yylex (void) {
  string s;
  int c;
  int j;
  
  c = getchar();
  // skip whitespace 
  while (c == ' ' || c == '\t') {c = getchar (); }

  if (isalpha (c)) { 
    s="";
    j=0; 
    while (isalpha(c) || isdigit(c) || c=='_') {
      s.push_back(c);
      c = getchar();
    }
    ungetc(c,stdin);
    yylval.string = new char[s.size()+1];
    strcpy(yylval.string,s.c_str());
    return IDENT;
  } else if (isdigit (c)) { 
    int len;
    s="";
    j=0; 
    while (isdigit(c)) {
      s.push_back(c);
      c = getchar();
    }
    ungetc(c,stdin);
    sscanf(s.c_str(),"%d",&len);
    yylval.num=len;
    return LENGTH;
  } else if (c == '[' || c == ']') {
    return c;
  } else if (c == '\n') {
    return 0;
  }
  return 1;
}

int main() {
  while (1) {
    prop_name.clear();
    prop_len.clear();
    printf("> ");
    yyparse();
    for (int j=0; j<prop_len.size(); j++) {
      printf("name %s, len %d\n",prop_name[j].c_str(),prop_len[j]);
    }
  }

}

/* Called by yyparse on error */
int yyerror (char *s) {
  printf ("%s\n", s);
  return 0;
}
