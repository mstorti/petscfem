/*__INSERT_LICENSE__*/
// $Id: tryme4.cpp,v 1.1 2002/07/18 23:38:54 mstorti Exp $

#include <cmath>
#include <deque>

double drand() {  
  return ((double)(rand()))/((double)(RAND_MAX));
}

int irand(int imin,int imax) {
  return int(rint(drand()*double(imax-imin+1)-0.5))+imin;
}

template<class T>
class SET {
private:
  typdef deque<T> cont;
  typdef cont::iterator cont_it;
  cont d;
  cont_it end_ord;
  // will resort if size passes max
  int max;
public:
  SET() { end_ord = d.end(); max=1000; }
  ~SET() { d.clear(); }
  insert(T t) {
    cont_it q = binary_search(d.begin(),end_ord,t);
    if (q==end_ord) {
      d.push_back(t);
      if (d.size()>max) {
	sort(d.begin(),d.end());
	max = 2*max;
      }
    }
  }
};

int main(int argc, char **argv) {
  SET<int> g;
}
