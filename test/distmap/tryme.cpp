//__INSERT_LICENSE__
//$Id: tryme.cpp,v 1.4 2003/01/08 15:49:03 mstorti Exp $
#include <iostream>
#include <list>

template<class B,class Cond>
class Grep {
public:
  Cond c;
  void apply(const list<B> &in, list<B> &out);
};
  

template<class B,class Cond>
void Grep<B,Cond>::apply(const list<B> &in, list<B> &out) {
  list<B>::const_iterator k;
  for (k=in.begin(); k!=in.end(); k++) {
    if (c.satisfies(*k)) {
      out.insert(out.end(),*k);
    }
  }
}

class Even {
public:
  int satisfies(int k) {return (k % 2 == 0);};
};

class Odd {
public:
  int satisfies(int k) {return (k % 2 != 0);};
};

int main() {
  Grep<int,Even> evengrep;
  Grep<int,Odd> oddgrep;
  list<int> in,out;
  list<int>::iterator j;

  for (int k=0; k<20; k++)
    in.insert(in.end(),k);

  evengrep.apply(in,out);

  cout << "Even grep: ";
  for (j=out.begin(); j!=out.end(); j++) 
    cout << *j << "  ";
  cout << endl;

  out.clear();
  oddgrep.apply(in,out);

  cout << "Odd grep: ";
  for (j=out.begin(); j!=out.end(); j++) 
    cout << *j << "  ";
  cout << endl;

}
