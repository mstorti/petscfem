// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: mmove.h,v 1.19 2005/06/11 13:11:56 mstorti Exp $

#ifndef MMOVE_H
#define MMOVE_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/** Distorts meshes solving a elasticity-like problem. 
    When the domain changes, you have to relocate nodes.
    At this time this can be done in PETSc-FEM only in a 
    separate run. In this run you enter the diplacements
    on the boundary as a fixation and using this elemenset 
    the resulting field is the displacement in the internal nodes. 
    Unlike the elasticity problems, this elemsets treats the mesh 
    relocations problem as the minimization of a functional of the
    distortion of the element w.r.t. a ``master'' (regular) element,
    (for instance an equilateral triangle, a square, a regular
    terahedra or a cube). We compute the metric tensor of the
    transformation from the actual element to the regular master element
    and the functional to be minimized is a function of the eigenvalues 
    of this tensor. 
 */
class  mesh_move_eig_anal : public adaptor { 
private:
  /// Auxiliary variables
  FastMat2 G, lambda, glambda, V, J, tmp1, tmp2, 
    dNdxi, xlocp, xloc0, x0,xp,lambdap,glp,glambda_diff,
    dFdl, d2Fdl2, resp, res_Dir, dstate, tmp3,tmp4;
  /// Parameters
  double eps,distor_exp,c_distor,c_volume,c_relax;
  /** Computes the eigenvalues of the meric tensor and its gradient. 
      @param x (input) the point (element nodes coordinates)
      @param lambda (output) the eigenvalues of the metric tensor at #x#
      @param glambda (output) the derivatives of the eigenvalues with respect
          to the element nodes coordinates. 
   */
  void la_grad(const FastMat2 &x,FastMat2 &lambda,
	       FastMat2 &glambda);
  /** Computes the gradient of the functional (residual) for
      a given node element coordinates. (This is done analytically). 
      @param x (input) node element coordinates
      @param dFdx (output) the derivatives (i.e. the gradient) of the
      functional w.r.t. the node element coordinates. 
   */ 
  void df_grad(const FastMat2 &x,FastMat2 &dFdx);
  /** The distortion function as a function of the eigenvalues. 
      @param D (input) the eigenvalues of the metric tensor. 
      @return the distortion functional (to be minimized). 
  */ 
  double dfun(const FastMat2 &D);
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
