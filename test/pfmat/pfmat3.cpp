/*__INSERT_LICENSE__*/
// $Id: pfmat3.cpp,v 1.4 2004/01/21 02:47:37 mstorti Exp $

#include <src/dvector.h>
#include <src/dvector2.h>
#include <cassert>
#include <cstdio>
#include <cmath>

#include "./hashf.h"

class pfmat2 {
private:
  class cell {
    int i,j;
    double val;
    friend class pfmat2;
    int next;
  };
  dvector<cell> cells;
  dvector<int> headers;
  int N;			// number of buckets
  int M;			// size of cell vector
  int top;			// first free cell
  int mask;
  int hash_fun(int i,int j) {
    // return ((i+j)*1324) % N;
    // return i % N;
    int h = ::hash_fun((ub1 *)&i,4,0);
    h = ::hash_fun((ub1 *)&j,4,h);
    printf("h(%d) = %d\n",j,h);
    return h & mask;
  };
  int free() {
    int c = top;
    assert(c>=0);
    top = cells.ref(c).next;
    return c;
  }
public:
  pfmat2(int N_a,int M_a=0) : N(N_a), M(M_a) {
    assert(N>0);
    if (!M) M=N;
    cells.mono(N);
    headers.mono(N);
    top = 0;
    for (int j=0; j<M-1; j++) 
      cells.ref(j).next = j+1;
    cells.ref(M-1).next = -1;
    for (int j=0; j<M; j++) 
      headers.ref(j) = -1;
    double v = log2(double(N));
    int nbits = int(v);
    assert(v == double(nbits));
    mask = hashmask(nbits);
  }
  void add(int i,int j, double val) {
    int bucket = hash_fun(i,j);
    int *cursor = &headers.ref(bucket);
    while (*cursor>=0) {
      cell & c = cells.ref(*cursor);
      if (c.i==i && c.j==j) break;
      cursor = &c.next;
    }
    if (*cursor<0) { 
      *cursor = free();
      headers.ref(bucket) = *cursor;
      cell & c = cells.ref(*cursor);
      c.i = i;
      c.j = j;
      c.val = val;
      c.next = -1;
    } else {
      cell & c = cells.ref(*cursor);
      c.val += val;
    }
  }
  void print() {
    for (int j=0; j<N; j++) {
      int cursor = headers.ref(j);
      while (cursor>=0) {
	cell & c = cells.ref(cursor);
	printf("%d, %d: %f\n",c.i,c.j,c.val);
	cursor = c.next;
      }
    }
  }
};

int main() {
  pfmat2 A(1024);
  for (int j=0; j<10; j++) 
    A.add(j,j,double(2*j));
  for (int j=0; j<10; j++) 
    A.add(j,j,double(2*j));
  A.print();
}

