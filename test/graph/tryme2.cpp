/*__INSERT_LICENSE__*/
// $Id: tryme2.cpp,v 1.3 2002/07/17 22:45:30 mstorti Exp $

#include <glib.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>

double drand() {  
  return ((double)(rand()))/((double)(RAND_MAX));
}

int irand(int imin,int imax) {
  return int(rint(drand()*double(imax-imin+1)-0.5))+imin;
}

int int_compare(const void *aa,const void *bb) {
  const int *a = (const int *)aa;
  const int *b = (const int *)bb;
  return (*a<*b) - (*a>*b);
}

int main() {
  GTree* tree = g_tree_new(int_compare);
  for (int j=0; j<100; j++) { 
    int k,l;
    k = irand(1,10);
    l = irand(1,10);
    printf("%d -> %d\n",k,l);
    g_tree_insert(tree,new int(k),new int(l)); 
  }
  for (int k=-10; k<20; k++) {
    int *v = (int *)g_tree_lookup(tree,&k);
    if (v) {
      printf("%d -> %d\n",k,*v);
    } else {
      printf("%d -> <UNDEFINED>\n",k);
    }
  }
}
