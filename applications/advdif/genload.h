// -*-mode: c++ -*-
//__INSERT_LICENSE__
//$Id: genload.h,v 1.6 2003/01/08 15:54:25 mstorti Exp $
#ifndef GENLOAD_H
#define GENLOAD_H

#define FASTMAT2SHELL FastMat2Shell_t
class FASTMAT2SHELL {
public:
  virtual void prod(FastMat2 &Ax, FastMat2 &x) {
    PETSCFEM_ERROR0("Not overloaded prod for this FastMat2Shell.");
  }
  void add(FastMat2 &S) {
    PETSCFEM_ERROR0("Not overloaded 'add' for this FastMat2Shell.");
  }
  virtual void init() {};
  //  virtual ~FASTMAT2SHELL()=0;
};

/// Generic surface flux function (film function) element
class LinearHFilmFun : public HFilmFun {
private:
  FastMat2 dU;

  class H;
  class HFull;
  class S;
  class SFull;
  class SNull;

  friend class H;
  friend class HFull;
  friend class S;
  friend class SFull;
  friend class SNull;

  class H : public FASTMAT2SHELL {
  public:
    LinearHFilmFun* l;
    virtual void prod(FastMat2 &Ax, FastMat2 &x)=0;
    virtual void jac(FastMat2 &A)=0;
    virtual void init() {};
    virtual void element_hook(ElementIterator &element) {};
    H(LinearHFilmFun *l_) : l(l_) {};
  };
  
  class HFull : public H {
  private:
    FastMat2 HH;
  public:
    void prod(FastMat2 &Ax, FastMat2 &x) {Ax.prod(HH,x,1,-1,-1);};
    void jac(FastMat2 &A);
    void init();
    void element_hook(ElementIterator &element);
    HFull(LinearHFilmFun *l) : H(l) {};
  };

  class HNull : public H {
  public:
    void prod(FastMat2 &Ax, FastMat2 &x) {Ax.set(0.);};
    void jac(FastMat2 &A) {A.set(0.);};
    HNull(LinearHFilmFun *l) : H(l) {};
  };

  class S : public FASTMAT2SHELL {
  public:
    LinearHFilmFun* l;
    virtual void add(FastMat2 &flux)=0;
    virtual void init() {};
    virtual void element_hook(ElementIterator &element) {};
    S(LinearHFilmFun *l_) : l(l_) {};
  };

  class SFull : public S {
  private:
    FastMat2 SS;
  public:
    void add(FastMat2 &flux) {flux.add(SS);};
    void init() {SS.resize(1,l->ndof);};
    void element_hook(ElementIterator &element);
    SFull(LinearHFilmFun *l_) : S(l_) {};
  };
  
  class SNull : public S {
  public:
    void add(FastMat2 &flux) {};
    SNull(LinearHFilmFun *l_) : S(l_) {};
  };
  
  int nel, ndof, nelprops;
  H *h;
  S *s;
  Property hfilm_coeff_prop, 
    hfilm_source_prop;
public:
  void q(FastMat2 &uin,FastMat2 &uout,FastMat2 &flux,
	 FastMat2 &jacin,FastMat2 &jacout);
  void q(FastMat2 &uin,FastMat2 &flux,FastMat2 &jacin);
  void init();
  void element_hook(ElementIterator &element);
  LinearHFilmFun(GenLoad *e) : HFilmFun(e) {};
  ~LinearHFilmFun();
};


/// Linear surface flux element
class lin_gen_load : public GenLoad { 
public: 
  LinearHFilmFun linear_h_film_fun;
  lin_gen_load() : linear_h_film_fun(this) {h_film_fun = &linear_h_film_fun;};
  ~lin_gen_load() {};
};

#endif
