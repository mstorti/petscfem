// -*- mode: C++ -*- 
//__INSERT_LICENSE__
//$Id: getarray.h,v 1.6 2001/04/15 22:39:03 mstorti Exp $
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
  void *add_block(void *,void *);
  void *create_list(void *ident_);
  void *f_id_subs_l(void *ident_,void *subsl);
  void *list_product(void *sl1,void *sl2);
  void *create_list_from_length(int length);
  void *dotify(void *sl);
#if 0
  int create_ident_l(char *ident,int length);
  int create_ident_list(char *ident,int list);
#endif
#ifdef __cplusplus
}
#endif
#endif
