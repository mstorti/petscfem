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
}

%token <string> IDENT
%token <num> LENGTH
     
%%
/* GRAMMAR RULES */

input:  /* empty */
      | input opt_def
;

opt_def: IDENT '[' LENGTH ']' { add_entry($1,$3);
                                /* printf("name %s, length %d\n",$1,$3); */
                                /* prop_table_add($1,$3); */
                              }
        | IDENT               { add_entry($1,1);
                                /* printf("name %s, length 1\n",$1); */
                               }
;

%%
/* ADDITIONAL C CODE */

