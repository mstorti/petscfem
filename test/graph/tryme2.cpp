/*__INSERT_LICENSE__*/
// $Id: tryme2.cpp,v 1.4 2002/07/18 02:57:33 mstorti Exp $

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

typedef void (*GRemoveFunc)(gpointer key,gpointer val);

void remove(gpointer key,gpointer val) {
  int *k = (int *)key;
  int *v = (int *)val;
  printf("elim: %d (%p) -> %d (%p)\n",*k,k,*v,v);
  delete k;
  delete v;
}

void g_tree_remove_key(GTree *tree,gpointer key,GRemoveFunc fun) {
  gpointer val = g_tree_lookup(tree,key);
  if (val) fun(key,val);
  g_tree_remove(tree,key);
}

int main() {
  GTree* tree = g_tree_new(int_compare);
  for (int j=0; j<100; j++) { 
    int k,l, *pk , *pl;
    k = irand(1,10);
    l = irand(1,10);
    pk = new int(k);
    pl = new int(l);
    printf("%d (%p) -> %d (%p)\n",k,pk,l,pl);
    g_tree_remove_key(tree,&k,&remove);
    g_tree_insert(tree,pk,pl); 
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
