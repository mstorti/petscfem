/*__INSERT_LICENSE__*/
// $Id: tryme4.cpp,v 1.2 2002/07/18 23:56:09 mstorti Exp $

#include <cstdio>
#include <cmath>
#include <deque>
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
    sort(d.begin(),d.end()); 
  }
public:
  SET() { end_ord = d.end(); max=1000; }
  ~SET() { d.clear(); }
  void insert(T t) {
    if (!binary_search(d.begin(),end_ord,t)) {
      d.push_back(t);
      if (d.size()>max) {
	sort_m();
	max = 2*max;
      }
    }
  }
  void print() { sort_m(); print2(); }
  void print2();
};

void SET<int>::print2() {
  for (int j=0; j<d.size(); j++) {
    printf("%d ",d[j]);
  }
  printf("\n");
}

int main(int argc, char **argv) {
  SET<int> g;
  for (int j=0; j<100; j++) { g.insert(irand(1,100)); }
  g.print();
}
