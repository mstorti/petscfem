// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: genload.h,v 1.2 2002/12/16 17:26:03 mstorti Exp $
#ifndef GENLOAD_H
#define GENLOAD_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Generic surface flux element
class GenLoad : public Elemset { 
public:
  /** These are to pass the state of the `H' quantities on the
      internal and external layers. `Hin' is an alias for `H'
  */
  FastMat2 H_m,H_out_m;
protected:
  int nel2;
public: 
  const FastMat2 &H,&H_out,&H_in;
  GenLoad() : H_in(H_m), H(H_m), H_out(H_out_m) {}
  virtual void start_chunk() {}
  virtual void element_hook(int element) {}
  // for one layer
  virtual void q(FastMat2 &uin,FastMat2 &flux,FastMat2 &jac);
  // for two layers
  virtual void q(FastMat2 &uin,FastMat2 &uout,
		 FastMat2 &flux_in, FastMat2 &flux_out,
		 FastMat2 &jac);
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Surface flux element where #flux_in = -flux_out#
class ConsGenLoad : public GenLoad { 
private:
  FastMat2 jac_aux;
public: 
  void start_chunk();
  virtual void start_chunk_c() {}
  void q(FastMat2 &uin,FastMat2 &uout,
	 FastMat2 &flux_in, FastMat2 &flux_out, FastMat2 &jac);
  // for two layers
  virtual void q(FastMat2 &uin,FastMat2 &uout,
		 FastMat2 &fluxin, FastMat2 &jac)=0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class lin_gen_load : public ConsGenLoad {
public:
  void q(FastMat2 &uin,FastMat2 &uout, FastMat2 &fluxin, FastMat2 &jac);
};

#endif
