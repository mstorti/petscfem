#include <iostream>
#include <list>

template<class B,class Cond>
class Grep {
public:
  Cond *c;
  void apply(const list<B> &in, list<B> &out);
  Grep(Cond *cc) : c(cc) {};  hash
};
  
template<class B,class Cond>
void Grep<B,Cond>::apply(const list<B> &in, list<B> &out) {
  list<B>::const_iterator k;
  for (k=in.begin(); k!=in.end(); k++) {
    if (c->satisfies(*k)) {
      out.insert(out.end(),*k);
    }
  }
}

template<>
Grep<int>::satisfies(int &k) {
  return c->satisfies(k);
}

class Even : public Cond<int> {
public:
  int satisfies(int k) {return (k % 2 == 0);};
};

int main() {
  Even e;
  Grep<int> evengrep(&e);
  list<int> in,out;
  list<int>::iterator j;

  for (int k=0; k<20; k++)
    in.insert(in.end(),k);

  evengrep.apply(in,out);

  for (j=out.begin(); j!=out.end(); j++) 
    cout << *j << "  ";

  cout << endl;

}
