/*__INSERT_LICENSE__*/
// $Id: tryme4.cpp,v 1.6 2002/07/19 13:44:59 mstorti Exp $

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

template<class T>
class SET {
private:
  typedef vector<T> cont;
  typedef cont::iterator cont_it;
  cont d;
  cont_it end_ord;
  // will resort if size passes max
  int max;
  // flag whether new elements have been inserted or not
  int modif;
  void resync() { 
    if (modif) {
      printf("before: ");
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
      modif = 0;
      printf("after: ");
      print2();
      printf("resyncing, size %d\n",d.size());
    }
  }
public:
  SET() { end_ord = d.end(); max=10; modif=0; }
  ~SET() { d.clear(); }
  void insert(T t) {
    if (!binary_search(d.begin(),end_ord,t)) {
      d.push_back(t); modif=1;
      if (d.size()>max) {
	resync();
	int new_max = 2*d.size();
	if (new_max>max) {
	  max = new_max;
	  printf("new size %d\n",max);
	}
      }
    } else printf("already in set...\n");
  }
  void print(const char *s=NULL) { 
    resync(); 
    if (s) printf("%s",s);
    print2(); 
  }
  void print2();
  int size() { resync(); return d.size(); }
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
  int M=10;
  int l=0;
  for (int j=0; j<2; j++) { 
    for (int kk=0; kk<500; kk++) { 
      l++;
      k = irand(1,1000);
      printf("insertando: %d\n",k);
      g.insert(k); 
      gg.insert(k); 
    }
    if (g.size()!=gg.size()) {
      printf("on l=%d bad: g.size(): %d, gg.size() %d\n",l, g.size(),gg.size());
      g.print("usando my_set<int>: ");
      print_set(gg,"usando set<int>: ");
      exit(0);
    } else printf("on l=%d OK: g.size(): %d, gg.size() %d\n",l, g.size(),gg.size());

  }

  g.print("usando my_set<int>: ");
  print_set(gg,"usando set<int>: ");

  printf("g.size(): %d, gg.size() %d\n",g.size(),gg.size());
}
