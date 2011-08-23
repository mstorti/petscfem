// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
// $Id: absolay.h,v 1.2 2006/04/18 02:16:56 mstorti Exp $
#ifndef PETSCFEM_ABSOLAY_H
#define PETSCFEM_ABSOLAY_H

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/elemset.h>
#include "./advective.h"
#include <src/hook.h>
#include <src/texthf.h>

class AbsorbingLayer;
extern AbsorbingLayer  *absorbing_layer_elemset_p;

/** Perfectly Matched Layer type elemset */
class AbsorbingLayer : public NewElemset { 
protected:
  /** A pointer to the flux function. Different
      physics are obtained by redefining the flux function
  */
  NewAdvDifFF *adv_diff_ff;

private:
  int flag;
  FastMat2 Habso,C,Ay,Hm,Uref;
  // store this amount of time steps
  int nstep_histo;
  int nnod,ndof;
  // U, adn W values are stored here
  dvector<double> uhist,whist,u,w;
  int Nx,Ny,nsaverotw;
  double hy,Dt;
  int use_addhoc_surface_abso;

public:
  /// Contructor from the pointer to the fux function
  AbsorbingLayer(NewAdvDifFF *adv_diff_ff_a=NULL) :
    adv_diff_ff(adv_diff_ff_a), flag(0) { } 
  ~AbsorbingLayer() {delete adv_diff_ff;}
  
  /// The assemble function for the elemset. 
  NewAssembleFunction new_assemble;
  /// The ask function for the elemset. 
  ASK_FUNCTION;

  void initialize();
  void init_hook();
  void time_step_pre(int step);
  void time_step_post(int step);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class abso_hook : public Hook {
private:
public:
  abso_hook() { }
  void init(Mesh &mesh_a,Dofmap &dofmap,const char *name) {
    assert(absorbing_layer_elemset_p);
    absorbing_layer_elemset_p->init_hook();
  }
  void time_step_pre(double time,int step) {
    absorbing_layer_elemset_p->time_step_pre(step);
  }
  void time_step_post(double time,int step,
		      const vector<double> &gather_values) {
    absorbing_layer_elemset_p->time_step_post(step);
  }
  void stage(const char *jobinfo,int stage,
	     double time,void *data) { }
  void close() { }
};

#endif
