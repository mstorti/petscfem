//__INSERT_LICENSE__
//$Id: parsefldm.cpp,v 1.1 2001/04/14 00:03:10 mstorti Exp $

#include <cstdio>
#include <ctype.h>
#include <cmath>

#include <string>
#include <vector>

#include "getarray.h"
#include "getarrgr.tab.h"
#include "parsefld.tab.h"

vector<int> *prop_len;
vector<string> *prop_name;
const char *line_to_parse;
vector<string> *node_data_field_list;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void add_entry()"
void add_entry(char *s,int l) {
  prop_len->push_back(l);
  prop_name->push_back(string(s));
  // printf("deleting %p with %s in add_entry\n",s,s);
  delete[] s;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int yylex ()"
int yylex (void) {
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
    yylval.string = new char[s.size()+1];
    strcpy(yylval.string,s.c_str());
    // printf("creating %p with %s\n",yylval.string,yylval.string);
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
    yylval.num=len;
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
#define __FUNC__ "void parse_props_line()"
void parse_props_line(const char *line,vector<string> &prop_name_,vector<int>
		      &prop_len_) {
  prop_name = &prop_name_;
  prop_len = &prop_len_;
  prop_name->clear();
  prop_len->clear();
  line_to_parse = line;
  yyparse();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void parse_fields_line(const char *line)"
void parse_fields_line(vector<string> &fields,const char *line) {
  node_data_field_list = &fields;
  line_to_parse = line;
  field_parse();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int field_lex (void)"
int field_lex (void) {
#if 0
  return yylex();
#endif
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
    yylval.string = new char[s.size()+1];
    strcpy(yylval.string,s.c_str());
    // printf("creating %p with %s\n",yylval.string,yylval.string);
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
    yylval.num=len;
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
#undef __FUNC__
#define __FUNC__ ""
/* Called by yyparse on error */
int yyerror (char *s) {
  printf ("parse props line error: %s\n", s);
  return 0;
}

// This is a test
#if 1
int main() {
  char input[1000];
  vector<string> fields;
  while (1) {
    printf("enter line to parse> ");
    scanf("%s",input);
    parse_fields_line(fields,input);
    for (int j=0; j<fields.size(); j++) {
      printf("%d -> %s\n",j+1,fields[j].c_str());
    }
  }
}
#endif
