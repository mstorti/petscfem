#include <cstdio>
#include <cmath>
#include <map>
#include "../../src/texthash.h"
#include "../../src/utils.h"

int main () {

  map<int,int> table;
  TextHashTable thash;
  thash.register_name(string("test table"));
  char k[20],v[20];
  int NN=100,nn=100,Nget=100;

  printf("Loading TextHashTable... %d elements in range 1,%d\n",NN,nn);
  for (int j=0; j<NN; j++) {
    int n=irand(1,nn);
    sprintf(k,"k%d",n);
    sprintf(v,"v%d",n);
    printf("adding %s -> %s\n",k,v);
    table[n]=0;
    thash.add_entry(k,v);
  }

  printf("Accessing elements %d times...\n",Nget);
  const char *vv;
  for (int j=0; j<Nget; j++) {
    int n=irand(1,nn);
    sprintf(k,"k%d",n);
    thash.get_entry(k,vv);
    if (vv) {
      assert(table.find(n) != table.end());
      table[n]++;
      // printf("acessing %s\n",vv);
    }
  }
  // thash.read(fstack);
  thash.print();
  map<int,int>::iterator kk;
  int OK=1,na;
  for (kk=table.begin(); kk!=table.end(); kk++) {
    int n=kk->first;
    printf("k%d accessed %d times\n",kk->first,kk->second);
    sprintf(k,"k%d",n);
    na = thash.access_count(k);
    if (na != kk->second) {
      OK=0;
      printf("counting for n=%d, not OK. Counted %d, in thash: %d\n",n,kk->second,na);
    }
  }
  printf(OK ? "Global counting OK\n" : "Global counting not OK\n");

  TextHashTable::print_stat();
}
