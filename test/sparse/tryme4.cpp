/*__INSERT_LICENSE__*/
// $Id: tryme4.cpp,v 1.9 2002/07/19 20:16:22 mstorti Exp $

#include <cassert>
#include <cstdio>
#include <cmath>
#include <deque>
#include <vector>
#include <set>
#include <algorithm>

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
class Set {
private:
  typedef deque<T,malloc_alloc> cont;
  typedef cont::iterator cont_it;
  cont d;
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
      sort(d.begin(),d.end());
      cont_it p=d.begin(), e=d.end(), q;
      if (p==e) {
	ordered = 0;
	return;
      }
      q=p;
      while (++q!=e) if (*q!=*p) *++p = *q;
      d.erase(++p,e);
      ordered = p-d.begin();
      modif = 0;
      // printf("after: ");
      // print2();
      // printf("resyncing, size %d\n",d.size());
    }
  }
public:
  Set() { ordered = 0; MAX_INIT = 1000; max=MAX_INIT; modif=0; }
  ~Set() { d.clear(); }
  void insert(T t) {
    if (!binary_search(d.begin(),d.begin()+ordered,t)) {
      d.push_back(t); modif=1;
      if (d.size()>max) {
	resync();
	int new_max = 2*d.size();
	if (new_max>max) {
	  max = new_max;
	  // printf("new size %d\n",max);
	}
      }
    } // else printf("already in set...\n");
  }
  void print(const char *s=NULL) { 
    resync(); 
    if (s) printf("%s",s);
    print2(); 
  }
  void print2();
  int size() { resync(); return d.size(); }
  void clear() { d.clear(); modif=0; ordered=0; max = MAX_INIT; }
};

void Set<int>::print2() {
  for (int j=0; j<d.size(); j++) {
    printf("%d ",d[j]);
  }
  printf("\n");
}

void Set<int_pair>::print2() {
  for (int j=0; j<d.size(); j++) {
    printf("(%d,%d) ",d[j].i,d[j].j);
  }
  printf("\n");
}

int main(int argc, char **argv) {
  Set<int_pair> g;
  Set2 gg;
  int kk;
  for (kk=0; kk<1000000; kk++) { 
    int k = irand(1,1000);
    int l = irand(1,1000);
    g.insert(int_pair(k,l)); 
    gg.insert(int_pair2(k,l)); 
    if (kk % 100000 == 0 ) {
      if (g.size()!=gg.size()) {
	printf("on kk=%d bad: g.size(): %d, gg.size() %d\n",kk, g.size(),gg.size());
	g.print("usando my_set<int>: ");
	print_set(gg,"usando set<int>: ");
	exit(0);
      } else printf("on kk=%d OK: g.size(): %d, gg.size() %d\n",kk, g.size(),gg.size());
    }
  }
  g.clear();
  gg.clear();
  printf("g.size(): %d, gg.size() %d\n",g.size(),gg.size());
}
