//__INSERT_LICENSE__
// $Id: trytree.cpp,v 1.4 2004/12/26 03:09:24 mstorti Exp $

#include <cstdlib>
#include <cassert>
#include <cmath>
#include <vector>
#include <iostream>
#include "./tree.h"

using namespace std;

int main1() {
  tree<int> A;
  typedef tree<int>::iterator node_t;
  node_t n = A.begin();
  n = A.insert(n,10);
  n = A.insert(n.lchild(),13);
  n = A.insert(n,12);
  node_t n12 = n;
  n = A.insert(n,11);
  n = A.insert(n12.lchild(),23);
  n = A.insert(n,22);
  n = A.insert(n,21);
  A.lisp_print();
  printf("\n");

  printf("Path from 23 to root:\n");
  n = A.find(23);
  while (n!=A.end()) {
    printf("%d ",*n);
    n = n.father();
  }
  printf("\n");
}

double drand() {
  return double(rand())/double(RAND_MAX);
}

void make_random_tree(tree<int> &T,tree<int>::iterator n,
		      int M,int level,int siblings) {
  double lambda = 1.0/(double(siblings)/double(level)+1.0);
  n=T.insert(n,rand() % M);
  tree<int>::iterator c=n.lchild();
  while (true) {
    if (drand()<lambda) break;
    make_random_tree(T,c,M,level+1,siblings);
    c=n.lchild();
  }
}

void make_random_tree(tree<int> &T,int M,int siblings) {
  make_random_tree(T,T.begin(),M,1,siblings);
}

void ord_pre(tree<int> &T,tree<int>::iterator n,
	     vector<int> &v) {
  if (n==T.end()) return;
  v.push_back(*n);
  tree<int>::iterator q = n.lchild();
  while (q!=T.end()) ord_pre(T,q++,v);
}

void ord_pre(tree<int> &T,vector<int> &v) {
  v.clear();
  ord_pre(T,T.begin(),v);
}

#define BIG 10000

void swtch(int &v) {
  v = (v<BIG ? v+BIG : v-BIG);
}

int main() {
  tree<int> A;
  make_random_tree(A,20,3);
  A.lisp_print();

  vector<int> op;
  ord_pre(A,op);
  printf("pre-order: ");
  for (int j=0; j<op.size(); j++) {
    printf("%d ",op[j]);
    tree<int>::iterator n = A.find(op[j]);
    assert(n!=A.end());
    swtch(*n);
  }
  printf("\n");

  for (int j=0; j<op.size(); j++) {
    tree<int>::iterator n = A.find(op[j]+BIG);
    assert(n!=A.end());
    swtch(*n);
  }

  
  for (int j=0; j<op.size(); j++) {
    printf("path to %d: ",op[j]);
    tree<int>::iterator n = A.find(op[j]);
    tree<int>::iterator q=n;
    assert(n!=A.end());
    
    while (q!=A.end()) {
      printf("%d ",(*q<BIG ? *q : *q-BIG));
      q = q.father();
    }
    printf("\n");

    swtch(*n);
    printf("tree %p, tree from cell %p\n",
	   n.tree_owner(),n.cell_ptr()->tree_owner());
  }

  tree<int> B;
  B.splice(B.begin(),A.begin());

  for (int j=0; j<op.size(); j++) {
    printf("path to %d: ",op[j]);
    tree<int>::iterator n = B.find(op[j]+BIG);
    assert(n != B.end());
    tree<int>::iterator q = n;
    
    while (q != B.end()) {
      printf("%d ",(*q < BIG ? *q : *q-BIG));
      q = q.father();
    }
    printf("\n");

    swtch(*n);
    printf("tree %p, tree from cell %p\n",
	   n.tree_owner(),n.cell_ptr()->tree_owner());
  }
}
