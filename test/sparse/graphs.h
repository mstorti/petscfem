// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: graphs.h,v 1.1 2002/07/21 16:00:36 mstorti Exp $
#ifndef GRAPHS_H
#define GRAPHS_H

double drand();
int irand(int imin,int imax);

class graph {
public:
  virtual void add(int j,int k)=0;
  virtual void clear()=0;
  virtual int size()=0;
  virtual void print(const char *s= NULL)=0;
};

typedef pair<int,int> int_pair2;

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

#endif
