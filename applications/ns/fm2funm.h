#ifndef FM2FUNM_H
#define FM2FUNM_H

class FastMat2_fund {
private:
  int m;
  FastMat2 V,D,tmp,lambda,flambda;
public:
  void init(FastMat2 &A);
  void apply(const FastMat2 &A,FastMat2 &fA);
  virtual void f(const FastMat2 &D,FastMat2 &fD)=0;
};

double FastMat2_funm_ff(double x,void *a);

class FastMat2_funm : public FastMat2_fund {
private:
  void f(const FastMat2 &D,FastMat2 &fD);
public:
  virtual double f(double)=0;
};

#endif
