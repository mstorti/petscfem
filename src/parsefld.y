/* C DECLARATIONS */

%{
#include <stdio.h>
#include <ctype.h>
#include "getarray.h"
#define YYDEBUG 1

%}

/* BISON DECLARATIONS */

%union{
  int num;
  char *string;
  void *gen_ptr;
}

%token <string> IDENT
%token <num> LENGTH
%type <num> field_block_list field_block input
     
%%
/* GRAMMAR RULES */

input:  /* empty */ { $$ = 0;}
  | field_block_list { $$ = $1;}
;

field_block_list: field_block { printf("first field block in list\n");
                                $$=$1;}
| field_block_list field_block {add_block($1,$2);}
;

field_block: IDENT {$$ = create_list($1);}
         | IDENT '[' LENGTH ']' { $$ = create_ident_l($1,$3);}
         | IDENT '[' field_block_list ']' {$$ = create_ident_list($1,$3);}
;

%%
