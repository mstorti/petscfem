//__INSERT_LICENSE__
//$Id: iterator.cpp,v 1.2 2001/04/01 01:35:07 mstorti Exp $

#include <stdio.h>
#include <vector>

//class mv_iterator;
class mv_const_iterator;

class myvector : public vector<int> {
public:
  friend class mv_iterator;
  mv_const_iterator begin(void) const;
  mv_const_iterator end(void) const;
};

class mv_const_iterator {
protected:
  int position;
  const myvector *mv;
public:
  mv_const_iterator(const int p=0,
		    const myvector *mv_=NULL) : position(p), mv(mv_) {};
  mv_const_iterator &operator++(void);
  int operator!=(mv_const_iterator other) const;
  const int &operator*(void) const;
};

#if 0
class mv_iterator : public mv_const_iterator {
private:
  myvector *mv;
public:
  mv_iterator(const int p=0, myvector *mv_=NULL) : position(p), mv(mv_) {};
  int &operator*(void) const;
};
#endif

mv_const_iterator myvector::begin(void) const {
  return mv_const_iterator(0,this);
}

mv_const_iterator myvector::end(void) const {
  return mv_const_iterator(size(),this);
}

mv_const_iterator &mv_const_iterator::operator++(void) {
  position++;
  return *this;
}

int mv_const_iterator::operator!=(mv_const_iterator other) const {
  return position!=other.position;
}

const int &mv_const_iterator::operator*(void) const {
  return (*mv)[position];
}

int main() {
  myvector mv;
  mv_const_iterator i;
  for (int k=0; k<10; k++) {
    mv.push_back(k);
  }

  int j=0;
  for (i=mv.begin(); i!=mv.end(); ++i) {
    printf("j: %d, mv[j]: %d\n",j,*i);
    j++;
  }
}
