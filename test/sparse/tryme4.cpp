/*__INSERT_LICENSE__*/
// $Id: tryme4.cpp,v 1.10 2002/07/20 13:50:17 mstorti Exp $

#include <cassert>
#include <cstdio>
#include <cmath>
#include <deque>
#include <vector>
#include <set>
extern "C" {
#include <src/libretto.h>
#include <libretto/darray.h>
}

double drand() {  
  return ((double)(rand()))/((double)(RAND_MAX));
}

int irand(int imin,int imax) {
  return int(rint(drand()*double(imax-imin+1)-0.5))+imin;
}

void print_set(set<int> &gg,const char *s=NULL) { 
    if (s) printf("%s",s);
    for (set<int>::iterator q=gg.begin(); q!=gg.end(); q++) 
    printf("%d ",*q);
}

typedef pair<int,int> int_pair2;
typedef set<int_pair2,less<int_pair2>,malloc_alloc> Set2;
void print_set(Set2 &gg,const char *s=NULL) { 
    if (s) printf("%s",s);
    for (set<int_pair2>::iterator q=gg.begin(); q!=gg.end(); q++) 
    printf("(%d,%d) ",q->first,q->second);
}

class int_pair { 
public: 
  int i,j; 
  int_pair(int ii,int jj) { i=ii; j=jj; }
  bool operator<(int_pair q) const {
    if (i!=q.i) return i<q.i;
    else return j<q.j;
  }
  bool operator==(int_pair q) const {
    return i==q.i && j==q.j;
  }
};

template<class T>
int compare(const void *aa,const void *bb, void * d) {
  T * a = (T *) aa;
  T * b = (T *) bb;
  return - (*a < *b) + (*a > *b);
}

template<class T>
class Set {
private:
  Darray *da;
  int ordered;
  int MAX_INIT;
  // will resort if size passes max
  int max;
  // flag whether new elements have been inserted or not
  int modif;
  void resync() { 
    if (modif) {
      // printf("before: ");
      // print2();
      da_sort(da,compare<T>,NULL);
      printf("antes length: %d\n",da_length(da));
      int p=0, e = da_length(da), q;
      if (p==e) {
	ordered = 0;
	return;
      }
      q=p;
      while (++q!=e) 
	if (compare<T>(at(q),at(p),NULL)) 
	  da_set (da,++p,da_ref(da,q));
      ordered = p+1;
      // da_resize(da,ordered);
      printf("despues length: %d\n",da_length(da));
      modif = 0;
      // printf("after: ");
      // print2();
      // printf("resyncing, size %d\n",d.size());
    }
  }
  T* at(int q) { return (T*) da_ref(da,q); }
public:
  Set() { 
    da = da_create (sizeof(int_pair)); 
    ordered = 0; MAX_INIT = 1000; 
    max=MAX_INIT; modif=0; 
  }
  ~Set() { da_destroy(da); }
  void insert(T t) {
    if (da_bsearch(da,&t,&compare<T>,NULL) < 0) {
      da_append(da,&t); modif=1;      
      int ds = da_length(da);
      if (ds > max) {
	resync();
	int new_max = 2*ds;
	if (new_max>max) {
	  max = new_max;
	  printf("new size %d\n",max);
	}
      }
    }  else printf("already in set...\n");
  }
  void print(const char *s=NULL) { 
    resync(); 
    if (s) printf("%s",s);
    print2(); 
  }
  void print2();
  int size() { resync(); return da_length(da); }
  void clear() { da_resize(da,0); modif=0; ordered=0; max = MAX_INIT; }
};

void Set<int>::print2() {
  int ds = da_length(da);
  for (int j=0; j<ds; j++) {
    printf("%d ",*at(j));
  }
  printf("\n");
}

void Set<int_pair>::print2() {
  int ds = da_length(da);
  for (int j=0; j<ds; j++) {
    printf("(%d,%d) ",at(j)->i,at(j)->j);
  }
  printf("\n");
}

int main(int argc, char **argv) {
  Set<int_pair> g;
  Set<int_pair> gg;
  // Set2 gg;
  int kk;
  int N = 100000, M = int(N/10), NN = int(sqrt(N/4));
  for (kk=1; kk<=N; kk++) { 
    int k = irand(1,NN);
    int l = irand(1,NN);
    // printf("inserting (%d,%d)\n",k,l);
    g.insert(int_pair(k,l)); 
    gg.insert(int_pair(k,l)); 
    if (kk % M == 0 ) {
      if (g.size()!=gg.size()) {
	printf("on kk=%d bad: g.size(): %d, gg.size() %d\n",kk, g.size(),gg.size());
	g.print("usando my_set<int>: ");
	gg.print("usando my_set<int>: ");
	// print_set(gg,"usando set<int>: ");
	exit(0);
      } else printf("on kk=%d OK: g.size(): %d, gg.size() %d\n",kk, g.size(),gg.size());
    }
  }
  g.clear();
  gg.clear();
  printf("g.size(): %d, gg.size() %d\n",g.size(),gg.size());
}
