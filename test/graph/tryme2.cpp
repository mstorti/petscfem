/*__INSERT_LICENSE__*/
// $Id: tryme2.cpp,v 1.2 2002/07/17 22:33:54 mstorti Exp $

#include <cstdio>
#include <glib.h>

int int_compare(const void *aa,const void *bb) {
  const int *a = (const int *)aa;
  const int *b = (const int *)bb;
  return (*a<*b) - (*a>*b);
}

int main() {
  GTree* tree = g_tree_new(int_compare);
  for (int j=0; j<10; j++) { 
    g_tree_insert(tree,new int(j),new int(j*j)); 
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

