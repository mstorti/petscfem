// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: graphs.h,v 1.2 2002/07/21 19:01:48 mstorti Exp $
#ifndef GRAPHS_H
#define GRAPHS_H

#include <algorithm>

double drand();
int irand(int imin,int imax);

class graph {
public:
  virtual void add(int j,int k)=0;
  virtual void clear()=0;
  virtual int size()=0;
  virtual void print(const char *s= NULL)=0;
  virtual void set_ngbrs(int v,vector<int> &ngbrs)=0;
};


class graph_stl : public graph {
private:
  typedef pair<int,int> int_pair2;
  typedef set<int_pair2> gset;
  gset s;
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
  virtual void set_ngbrs(int v,vector<int> &ngbrs) {
    gset::iterator b,e,q;
    b = lower_bound(s.begin(), s.end(), int_pair2(v,0));
    e = lower_bound(s.begin(), s.end(), int_pair2(v+1,0));
    for (q=b; q!=e; q++) ngbrs.push_back(q->second);
  }
};

#endif
