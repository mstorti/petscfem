/*__INSERT_LICENSE__*/
//$Id: tautostr.cpp,v 1.1 2003/02/16 15:15:53 mstorti Exp $
#define _GNU_SOURCE
#include <cstdio>
#include <src/autostr.h>

int main () {
  AutoString s,t;
  int k=4,l=5,m=45,n=456;
  s.sprintf("<asjhas ajhs ajsh %d ajsh %d %d %d>",k,l,m,n);
  t.sprintf("<asjhas ajhs ajsh %d ajsh %d wewe we %d wewe wwewe we %d>",k,l,m,n);
  s.cat(t);
  s.cat_sprintf("<asjhas %d %d ajsh %d we %d>",k,l,m,n);
  printf("\"%s\"\n",s.str());
}
