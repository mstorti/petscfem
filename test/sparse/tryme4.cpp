/*__INSERT_LICENSE__*/
// $Id: tryme4.cpp,v 1.3 2002/07/19 00:31:39 mstorti Exp $

#include <cassert>
#include <cstdio>
#include <cmath>
#include <deque>
#include <set>
#include <algorithm>

double drand() {  
  return ((double)(rand()))/((double)(RAND_MAX));
}

int irand(int imin,int imax) {
  return int(rint(drand()*double(imax-imin+1)-0.5))+imin;
}

template<class T>
class SET {
private:
  typedef deque<T> cont;
  typedef cont::iterator cont_it;
  cont d;
  cont_it end_ord;
  // will resort if size passes max
  int max;
  void sort_m() { 
    print2();
    sort(d.begin(),d.end());
    cont_it p=d.begin(), e=d.end(), q;
    if (p==e) {
      end_ord = e;
      return;
    }
    q=p;
    while (++q!=e) if (*q!=*p) *++p = *q;
    d.erase(++p,e);
    end_ord = p;
    print2();
  }
public:
  SET() { end_ord = d.end(); max=10; }
  ~SET() { d.clear(); }
  void insert(T t) {
    if (!binary_search(d.begin(),end_ord,t)) {
      d.push_back(t);
      if (d.size()>max) {
	sort_m();
	max = 2*max;
	printf("resorting, max= %d\n",max);
      }
    }
  }
  void print() { sort_m(); print2(); }
  void print2();
  int size() { return d.size(); }
};

void SET<int>::print2() {
  for (int j=0; j<d.size(); j++) {
    printf("%d ",d[j]);
  }
  printf("\n");
}

int main(int argc, char **argv) {
  SET<int> g;
  set<int> gg;
  int k;
  for (int j=0; j<1000; j++) { 
    k = irand(1,1000);
    g.insert(k); 
    gg.insert(k); 
  }
  g.print();
  assert(g.size()==gg.size());
}
