// -*- mode: C++ -*- 
//__INSERT_LICENSE__
//$Id: getarray.h,v 1.5 2001/04/15 00:58:44 mstorti Exp $
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

  int field_parse(void);
  int field_lex(void);
  int field_error(char *s);
  void add_block(int,int);
  int create_list(void *ident_);
  int f_id_subs_l(void *ident_,int subsl);
  int list_product(int sl1,int sl2);
  int create_list_from_length(int length);
  int dotify(int sl);
#if 0
  int create_ident_l(char *ident,int length);
  int create_ident_list(char *ident,int list);
#endif
#ifdef __cplusplus
}
#endif
#endif
