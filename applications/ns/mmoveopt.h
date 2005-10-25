// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: mmoveopt.h,v 1.1 2005/10/25 12:43:53 mstorti Exp $

#ifndef MMOVEOPT_H
#define MMOVEOPT_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/** Distorts meshes solving a elasticity-like problem. 
    When the domain changes, you have to relocate nodes.
    At this time this can be done in PETSc-FEM only in a 
    separate run. In this run you enter the diplacements
    on the boundary as a fixation and using this elemenset 
    the resulting field is the disaplcement in the internal nodes. 
    Unlike the elasticity problems, this elemsets treats the mesh 
    relocations problem as the minimization of a functional of the
    distortion of the element w.r.t. a ``master'' (regular) element,
    (for instance an equilateral triangle, a square, a regular
    terahedra or a cube). We compute the metric tensor of the
    transformation from the actual element to the regular master element
    and the functional to be minimized is a function of the eigenvalues 
    of this tensor. 
 */
class  mesh_move_opt : public adaptor { 
private:
  /// Auxiliary variables
  FastMat2 dVdW,dSldW,dWdu,d2VdW2,d2SldW2,d2Vdu2,d2Sldu2,x,w,dVdu,dSldu,dQ,d2Q,tmp,mat1,vaux,vaux1,vaux2,w0,x0,epsilon_LC;
  /// Parameters
  double distor_exp,c_distor,c_volume,c_relax,volume_exp;

public: 
  /** Initializes the elemset. Reads parameters, 
      resizes auxiliary matrices. 
  */ 
  void init();
  /** Computes the residual (gradient of the distortion functional) and 
      the jacobian (Hessian, i.e. second derivatives, of the 
      distortion functional. 
      @param xloc (input) element node coordinates
      @param state_old (input) element node displacements (previous step)
      @param state_new (input) element node displacements (actual step)
      @param res (output) gradient of distortion functional
      @param mat (input) Hessian (matrix of second derivatives) 
      of distortion functional 
   */ 
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
};

#endif
