/* C DECLARATIONS */

%{
#define YYSTYPE char *
#include <stdio.h>
#include <ctype.h>
#define MAXLEN 20
%}

/* BISON DECLARATIONS */

%token IDENT
%token LENGTH
%token SPACE
     
%%
/* GRAMMAR RULES */

input:  /* empty */
      | opt_def
      | input SPACE opt_def
;

opt_def: IDENT '[' LENGTH ']' { printf("name %s, length %s\n",$1,$3);}
        | IDENT               { printf("name %s, length 1\n",$1);}
;

%%
	/* ADDITIONAL C CODE */

yylex ()
{
  int c;
  int j;
  
  c = getchar();
  if (c == ' ') 
    { 
      while (c == ' ') {c = getchar (); }
      ungetc(c,stdin);
      return SPACE;
    }
  else if (isalpha (c))
    { 
      yylval = (char *) malloc (MAXLEN * sizeof (char *));
      j=0; 
      while (isalpha(c) || isdigit(c) || c=='_') {
	yylval[j++] = c;
	c = getchar();
      }
      ungetc(c,stdin);
      yylval[j]='\0';
      return IDENT;
    }
  else if (isdigit (c))
    { 
      yylval = (char *) malloc (MAXLEN * sizeof (char *));
      j=0; 
      while (isdigit(c)) {
	yylval[j++] = c;
	c = getchar();
      }
      ungetc(c,stdin);
      yylval[j]='\0';
      return LENGTH;
    }
  else if (c == '[' || c == ']') 
    {
      return c;
    }
  else if (c == '\n') {
    return 0;
  }
}

int main() {
  while (1) {
    printf("> ");
    yyparse();
  }
}

yyerror (s)  /* Called by yyparse on error */
     char *s;
{
  printf ("%s\n", s);
}
