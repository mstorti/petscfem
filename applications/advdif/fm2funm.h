//__INSERT_LICENSE__
//$Id: fm2funm.h,v 1.1.4.1 2004/05/21 23:46:28 mstorti Exp $
#ifndef FM2FUNM_H
#define FM2FUNM_H

class FastMat2_fund {
private:
  int m;
  FastMat2 V,D,Vi,tmp,lambda,flambda;
public:
  void init(FastMat2 &A);
  void apply(const FastMat2 &A,FastMat2 &fA);
  void apply_ns(const FastMat2 &A,FastMat2 &fA);
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
