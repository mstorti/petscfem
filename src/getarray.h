// -*- mode: C++ -*- 
#ifndef GETARRAY_H
#define GETARRAY_H

#define MAXLEN 20

#ifdef __cplusplus
extern "C" {
#endif
void add_entry(char *s,int l);
int yyparse(void);
int yylex(void);
int yyerror(char *s);
#if 0
void add_index(void **sl,char *s,int last);
void opt_def_sl_action(char *s,void *sl);
#endif
#ifdef __cplusplus
}
#endif
#endif
