/*__INSERT_LICENSE__*/
// $Id: pfmat3.cpp,v 1.2 2004/01/20 22:44:11 mstorti Exp $

// Tests for the `PFMat' class

#if 1
#include <ext/hash_map>
#include <cstdio>

using namespace __gnu_cxx;
using namespace std;

int main(int argc,char **args) {
  hash<int> H;
  for (int j=0; j<1000; j++) {
    int h = H(j);
    printf("%d -> %d\n",j,h);
  }
}

#else
#include <glib.h>
#include <cstdio>
#include <cstring>

int my_int_hash_fun(int j) {
  char a[5];
  memcpy(a,&j,sizeof(int));
  a[4]='\0';
  return g_str_hash(a);
}

int main(int argc,char **args) {
  for (int j=0; j<1000; j++) {
    int h = my_int_hash_fun(j);
    printf("%d -> %d\n",j,h);
  }
}
#endif
