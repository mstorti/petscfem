/* C DECLARATIONS */

%{
#include <stdio.h>
#include <ctype.h>
#include "getarray.h"
%}

/* BISON DECLARATIONS */

%union{
  int num;
  char *string;
  void *gen_ptr;
}

%token <string> IDENT
%token <num> LENGTH
/* %type <gen_ptr> subscript_list */
     
%%
/* GRAMMAR RULES */

input:  /* empty */
      | input opt_def
;

opt_def: IDENT '[' LENGTH ']' { add_entry($1,$3);
                                /* printf("name %s, length %d\n",$1,$3); */
                                /* prop_table_add($1,$3); */
                              }
        | IDENT               { add_entry($1,1);}
/*| IDENT '[' subscript_list ']' { opt_def_sl_action($1,$3); } */
;

/* 
subscript_list: IDENT                   { add_index(&($$),$1,1);}
                | subscript_list IDENT  { add_index(&($$),$2,0);} 
;
*/

%%
/* ADDITIONAL C CODE */
