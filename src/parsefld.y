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
  void *gen_ptr;
}

%token <gen_ptr> IDENT
%token <num> LENGTH
%type <num> field_block_list field_block input subs_l subs
     
%%
/* GRAMMAR RULES */

input:  /* empty */ { $$ = 0;}
  | field_block_list { $$ = $1;}
;

field_block_list: field_block { $$=$1;}
       | field_block_list field_block {add_block($1,$2);}
;

field_block: IDENT {$$ = create_list($1);}
           | IDENT subs_l {$$ = f_id_subs_l($1,$2);}
;

subs_l: subs               { $$ = $1;}
	   | subs_l subs { $$ = list_product($1,$2);}
;

subs: '[' LENGTH ']'             {$$ = create_list_from_length($2);}
   |  '[' field_block_list ']' {$$ = dotify($2);}
;
  

%%
