// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: genload.h,v 1.6 2006/02/02 19:26:37 mstorti Exp $
#ifndef GENLOAD_H
#define GENLOAD_H

#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Generic surface flux element. See documentation 
    on the similar `advdif' elemset. In a future they
    will be unified. */
class GenLoad : public Elemset { 
public:
  /** These are to pass the state of the `H' quantities on the
      internal and external layers. `Hin' is an alias for `H' */
  FastMat2 H_m,H_out_m;
protected:
  int nel2;
  // Physical properties
  int elprpsindx[MAXPROP]; 
  int nprops;
  double propel[MAXPROP];
public: 
  const FastMat2 &H,&H_out,&H_in;
  GenLoad() : H_in(H_m), H(H_m), H_out(H_out_m) {}
  /** Call back function to be called by the elemset before
      a sequence of elements to be computed. This may be used
      by a derived class in order to perform some
      initializations. */
  virtual void start_chunk() {}
  /** Call back function to be called by the elemset after
      a sequence of elements to be computed. This may be used
      by a derived class in order to perform some
      cleanup. */
  virtual void end_chunk() {}
  /** fixme:= agregar doc 
      @param (input)
      @return a reference to the matrix. */ 
  virtual void element_hook(int element) {}
  /** One layer callback flux function. 
      @param u (input) state at the surface (size #ndof#)
      @param flux (output) flux to the surface (size #ndof#)
      @param jac (output) Jacobian of the flux with respect to the 
      vector state (#jac = (d flux)/(d u)#, size #ndox x ndof#) */ 
  virtual void q(FastMat2 &u,FastMat2 &flux,FastMat2 &jac);
  /** Two layer callback flux function. 
      @param uin (input) state at the internal surface (size #ndof#)
      @param uout (input) state at the exterior surface (size #ndof#)
      @param flux_in (output) flux from the exterior (size #ndof#)
      surface to the interior surface (size #ndof#)
      @param flux_out (output) flux from the interior
      surface to the exterior surface (size #ndof#)
      @param jac (output) Jacobian of the flux with respect to 
      the state of both surfaces (size #2 x 2 x ndof x ndof#, 
      where #jac(1,1,.,.)# is the jacobian #(d flux_in)/(d uin)#, 
      #jac(1,2,.,.) = (d flux_in)/(d u_out)#, 
      #jac(2,1,.,.) = (d flux_out)/(d u_in)#, 
      #jac(2,2,.,.) = (d flux_out)/(d u_out)#. */ 
  virtual void q(FastMat2 &u_in,FastMat2 &u_out,
		 FastMat2 &flux_in, FastMat2 &flux_out,
		 FastMat2 &jac);
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Surface flux element where #flux_in = -flux_out#. It's called
    conservative since the sum of the fluxes is null. (Normally this
    is true).  */
class ConsGenLoad : public GenLoad { 
private:
    /** Auxiliary variable, stores part of the jacobian to be copied */ 
  FastMat2 jac_aux;
  /// Initialization function.
  void start_chunk();
  /// Two layer callback flux function, implemented by this class. 
  void q(FastMat2 &uin,FastMat2 &uout,
	 FastMat2 &flux_in, FastMat2 &flux_out, FastMat2 &jac);
public: 
  /** Call back function to be called by the elemset before
      a sequence of elements to be computed. This may be used
      by a derived class in order to perform some
      initializations. */
  virtual void start_chunk_c() {}
  /** Two layer callback flux function. 
      @param uin (input) state at the internal surface (size #ndof#)
      @param uout (input) state at the exterior surface (size #ndof#)
      @param flux_in (output) flux from the exterior (size #ndof#)
      surface to the interior surface (size #ndof#)
      @param jac (output) Jacobian of the flux with respect to 
      the state of both surfaces (size #2 x ndof x ndof#, 
      where #jac(1,.,.)# is the jacobian #(d flux_in)/(d uin)#, 
      #jac(2,.,.) = (d flux_in)/(d u_out)#. */ 
  virtual void q(FastMat2 &uin,FastMat2 &uout,
		 FastMat2 &fluxin, FastMat2 &jac)=0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Flux function for a linear relation between the flux and the 
    surface states. */
class lin_gen_load : public ConsGenLoad {
private:
  ///  #h_film# Stores the matrix of coefficients: #flux = h_film (U_out-U_in)# 
  FastMat2 h_film, const_flux, U_out_sl;
  /// Temporary variable
  FastMat2 tmp1;
  int const_flux_indx;
  
public:
  /** Two layer callback flux function. 
      @param u_in (input) state at the internal surface (size #ndof#)
      @param u_out (input) state at the exterior surface (size #ndof#)
      @param flux_in (output) flux from the exterior (size #ndof#)
      surface to the interior surface (size #ndof#)
      @param Jac (output) is #Jac = [h_film -h_film]# */ 
  void q(FastMat2 &u_in,FastMat2 &u_out, 
	 FastMat2 &flux_in, FastMat2 &Jac);
  /// One layer callback flux function. 
  void q(FastMat2 &u,FastMat2 &flux,FastMat2 &jac);
  /** Initializes */ 
  void start_chunk_c();
};

#endif
