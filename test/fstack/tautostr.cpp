/*__INSERT_LICENSE__*/
//$Id: tautostr.cpp,v 1.4 2003/05/04 16:43:50 mstorti Exp $
#define _GNU_SOURCE
#include <cstdio>
#include <src/autostr.h>

int main () {
  AutoString s,t,u;
  int k=4,l=5,m=45,n=456,iter=0,niter=1;
  while (1) {
    s.sprintf("<This is string 1. A number: %d, other nubers: %d, %d, %d.>",k,l,m,n);
    t.sprintf("<This is string 2. Numbers: %d, %d, %d, %d.>",k,l,m,n);
    s.cat("\n")
      .cat(t)
      .cat("\n")
      .cat_sprintf("<String 3. Numbers %d, %d. Other %d, and %d.>",k,l,m,n)
      .cat("\n")
      .print();
    
    printf("\n\ns 2 times:\n");
    u.set(s)
      .cat("\n")
      .set(s)
      .cat("\n")
      .print();
    if (iter++==niter) break;
  }
}
