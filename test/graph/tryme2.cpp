/*__INSERT_LICENSE__*/
// $Id: tryme2.cpp,v 1.1 2002/07/17 22:15:25 mstorti Exp $

#include <glib.h>

int int_compare(const int *a,const in *b) {
  return a<b - a>b;
}

int main() {
  GTree* tree = g_tree_new(int_compare);
}

