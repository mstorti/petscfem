// -*- mode: C++ -*- 
//__INSERT_LICENSE__
//$Id: getarray.h,v 1.3 2001/04/01 01:35:03 mstorti Exp $
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
#ifdef __cplusplus
}
#endif
#endif
