// -*- mode: C++ -*- 
#ifndef GETARRAY_H
#define GETARRAY_H

#define MAXLEN 20

extern "C" {
void add_entry(char *s,int l);
int yyparse(void);
int yylex(void);
int yyerror(char *s);
}
#endif
