/*__INSERT_LICENSE__*/
// $Id: tryme4.cpp,v 1.15 2002/07/21 03:23:44 mstorti Exp $

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

template<class T>
class dvector {
private:
  // size of each chunk
  int chunk_size;
  // total number of elements 
  int size;
  // Number of chunks position of first free element
  int nchunks;
  // pointers to chunks
  vector<T *> chunk_vector;
  // pointer to STL storage array in chunk_vector
  T **chunks;
  void reff(int j,int &chunk,int &k) {
    chunk = j/chunk_size;
    k = j % chunk_size;
  }
public:
  dvector() { 
    chunk_size = 100;
    chunks = NULL;
    size = nchunks = 0;
  }
  T &ref(int j) {
    assert(j>=0 && j<size);
    int chunk,k;
    reff(j,chunk,k);
    return chunks[chunk][k];
  }
  void push(const T &t) {
    int chunk,k;
    reff(size,chunk,k);
    if (k==0) {
      assert(nchunks == chunk_vector.size());
      chunk_vector.push_back(new T[chunk_size]);
      chunks = chunk_vector.begin();
      nchunks++;
    }
    ref(size++) = t;
  }
};

class graph {
public:
  virtual void add(int j,int k)=0;
  virtual void clear()=0;
  virtual int size()=0;
  virtual void print(const char *s= NULL)=0;
};

typedef pair<int,int> int_pair2;

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
      // printf("antes length: %d\n",da_length(da));
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
      Darray *new_da = da_create_len (sizeof(T),ordered); 
      for (int j=0; j<ordered; j++) 
	*(T *)da_ref(new_da,j) = *(T *)da_ref(da,j);
      da_destroy(da);
      da = new_da;
      // printf("despues length: %d\n",da_length(da));
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
    }  // else printf("already in set...\n");
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

void Set<int_pair>::print2() {
  int ds = da_length(da);
  for (int j=0; j<ds; j++) {
    printf("(%d,%d) ",at(j)->i,at(j)->j);
  }
  printf("\n");
}

class graph_da : public graph {
private:
  Set<int_pair> sip;
public:
  void add(int j,int k) { sip.insert(int_pair(j,k)); }
  void clear() { sip.clear(); }
  void print(const char *s = NULL) { sip.print(s); }
  int size() { return sip.size(); } 
};  

class graph_stl : public graph {
private:
  set<int_pair2> s;
public:
  void add(int j,int k) { s.insert(int_pair2(j,k)); }
  void clear() { s.clear(); }
  int size() { return s.size(); } 
  void print(const char *ss=NULL) { 
    if (ss) printf("%s",ss);
    for (set<int_pair2>::iterator q=s.begin(); q!=s.end(); q++) 
      printf("(%d,%d) ",q->first,q->second);
    printf("\n");
  }
};

int main(int argc, char **argv) {
  dvector<int> v;
  int M=1000;
  for (int j=0; j<M; j++) v.push(j*j);
  for (int j=0; j<M; j++) {
    if (j % (M/10) ==0 ) 
      printf("j %d, j^2 %d, v(j) %d\n",j,j*j,v.ref(j));
    assert(j*j==v.ref(j));
  }

#if 0
  graph_da g;
  graph_stl gg;
  int kk;
  int N = 100000, M = int(N/10), NN = int(sqrt(N/2));
  for (kk=1; kk<=N; kk++) { 
    int k = irand(1,NN);
    int l = irand(1,NN);
    // printf("inserting (%d,%d)\n",k,l);
    g.add(k,l); 
    gg.add(k,l); 
    if (kk % M == 0 ) {
      if (g.size()!=gg.size()) {
	printf("on kk=%d bad: g.size(): %d, gg.size() %d\n",kk, g.size(),gg.size());
	g.print("g: ");
	gg.print("gg: ");
	exit(0);
      } else printf("on kk=%d OK: g.size(): %d, gg.size() %d\n",kk, g.size(),gg.size());
    }
  }
  g.clear();
  gg.clear();
  printf("g.size(): %d, gg.size() %d\n",g.size(),gg.size());
#endif
}
