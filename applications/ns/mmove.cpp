//__INSERT_LICENSE__
//$Id: mmove.cpp,v 1.17 2002/12/12 19:35:19 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/texthf.h>
#ifdef USE_NEWMAT
#include <newmatap.h>
#endif

#include "nsi_tet.h"
#include "adaptor.h"
#include "mmove.h"

extern GlobParam *GLOB_PARAM;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mesh_move_eig_anal::init() {
  int ierr;
  // The constitutive eq. is stress = C strain
  FastMat2 C;
  double c1 = sqrt(1./3.), c2 = sqrt(1./6.);
  // These are the gradient of shape functions for a master
  // tetra, with unit edge length with nodes at 
  // [+-1/2 0 0], [0,sqrt(3)/2,0], [0,1/sqrt(6),sqrt(2/3)]
  double c[12] = {-1., -c1, -c2, +1., -c1, -c2, 
		  0., 2*c1, -c2, 0., 0., 3*c2};
  
  G.resize(2,ndim,ndim); // the metric tensor
  xlocp.resize(1,nel*ndim); // The perturbed coordinates
  xloc0.resize(1,nel*ndim); // The perturbed coordinates
  assert(ndof==ndim);
  assert(ndim==2 || ndim==3);
  // assert(nel==ndim+1); // Only for triangles in 2D, tetras in 3D
  J.resize(2,ndim,ndim);
  dNdxi.resize(2,ndim,nel);

  if (ndim==2) {
    if (nel==3) {
      dNdxi.setel(-sin(M_PI/3)*cos(M_PI/6),1,1);
      dNdxi.setel(-sin(M_PI/3)*sin(M_PI/6),2,1);
      dNdxi.setel(+sin(M_PI/3)*cos(M_PI/6),1,2);
      dNdxi.setel(-sin(M_PI/3)*sin(M_PI/6),2,2);
      dNdxi.setel(0                       ,1,3);
      dNdxi.setel(+sin(M_PI/3)            ,2,3);
    } else if (nel==4) {
      double cq[] = {-1,-1,1,-1,1,1,-1,1};
      C.resize(2,nel,ndim).set(cq).t();
      dNdxi.set(C);
    } else PETSCFEM_ERROR("Only tringles ad quads in 2D: nel %d\n",nel);
  } else {
    C.resize(2,nel,ndim).set(c).t();
    dNdxi.set(C);
  }
  // gradient of eigenvalues
  glambda.resize(3,ndim,nel,ndim);
  // Gradient and Hessian of functional (w.r.t. eigenvalues)
  dFdl.resize(1,ndim);
  d2Fdl2.resize(2,ndim,ndim);

  //o Perturbation scale length for increment in computing
  // the Jacobian with finite differences. 
  TGETOPTDEF(thash,double,epsilon_x,1.e-4);
  eps = epsilon_x;
  //o The functional to be minimized is $\Phi = \sum_{e=1,...,Nel} \phi_e^r#,
  // where #\phi_e = \sum_{i\neq j} (\lambda_i-\lambda_j)^2/Vol^{2/n_d},
  // and $r={\tt distor\_exp}$. 
  TGETOPTDEF_ND(thash,double,distor_exp,1.);
  //o Adds a term $\propto {\tt c\_volume}\,{\rm volume}$ to the functionala. 
  TGETOPTDEF_ND(thash,double,c_volume,0.);
  //o Compute an initial ``predictor'' step with this relaxation scale. 
  TGETOPTDEF_ND(thash,double,c_relax,1.);
  //o Scales distortion function
  TGETOPTDEF_ND(thash,double,c_distor,1.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mesh_move_eig_anal::la_grad(const FastMat2 &x,FastMat2 &lambda,
				 FastMat2 &glambda) {
  // Computes jacobian of master element
  //        -> actual element transformation
  J.prod(dNdxi,x,2,-1,-1,1);
  double detJaco;
  detJaco = J.det();
  if (detJaco <= 0.) {
    PETSCFEM_ERROR("Jacobian of element %d is negative or null\n"
		   " Jacobian: %f\n",elem,detJaco);  
  }
  G.prod(J,J,-1,1,-1,2);
  lambda.seig(G,V);
  tmp3.prod(G,V,1,-1,-1,2);
  tmp4.prod(V,tmp3,-1,1,-1,2);
#define SHF(n) n.print(#n ": ")
#ifdef DEBUG_ANAL
  tmp4.print("V' G V (D?): ");
  SHF(J);
  SHF(G);
  SHF(lambda);
  SHF(V);
#endif
  tmp1.prod(J,V,1,-1,-1,2);
  tmp2.prod(dNdxi,V,-1,1,-1,2);
  for (int q=1; q<=ndim; q++) {
    glambda.ir(1,q);
    tmp1.ir(2,q);
    tmp2.ir(2,q);
    glambda.prod(tmp1,tmp2,2,1);
  }
  glambda.rs().scale(2.);
  tmp1.rs();
  tmp2.rs();
}

class Prod : public FastMat2::Fun2 {
public:
  void pre() { set(1.); }
  double fun2(double a,double val) { return val*a; }
} prod;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double mesh_move_eig_anal::dfun(const FastMat2 &D) {
#if 1
  double F=0;
  double vol=1.;
  for (int k=1; k<=ndim; k++) vol *= D.get(k);
  for (int k=2; k<=ndim; k++)
    for (int l=1; l<k; l++) F += square(D.get(k)-D.get(l));
  F /= pow(vol,2./double(ndim));
  return pow(F,distor_exp);
#elif 0
  double p=distor_exp;
  double norm_D = D.norm_p_all(p);
  double norm_iD = D.norm_p_all(-p);
  return norm_D*norm_iD;
#else
  return D.sum_all()*sqrt(D.assoc_all(prod));
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mesh_move_eig_anal::df_grad(const FastMat2 &x,FastMat2 &dFdx) {
  la_grad(x,lambda,glambda);
  double F, F0 = dfun(lambda);
  for (int k=1; k<=ndim; k++) {
    lambdap.set(lambda).addel(eps,k);
    F  = dfun(lambdap);
    dFdl.setel((F-F0)/eps,k);
  }
  dFdx.prod(dFdl,glambda,-1,-1,1,2);
#ifdef DEBUG_ANAL
    x.print("x:");
    lambda.print("lambda:");
    glambda.print("glambda:");
    dFdl.print("dFdl:");
    dFdx.print("dFdx:");
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mesh_move_eig_anal::
element_connector(const FastMat2 &xloc,
		  const FastMat2 &state_old,
		  const FastMat2 &state_new,
		  FastMat2 &res,FastMat2 &mat) {

  x0.set(xloc).add(state_old);
  df_grad(x0,res);
  x0.reshape(1,nel*ndim);
  mat.reshape(3,nel,ndim,nel*ndim);
  for (int k=1; k<=nel*ndim; k++) {
    xp.set(x0).addel(eps,k).reshape(2,nel,ndim);
    df_grad(xp,resp);
    xp.reshape(1,nel*ndim);
#ifdef DEBUG_ANAL
    printf("k %d\n",k);
    x0.print("x0:");
    xp.print("xp:");
    resp.print("resp:");
    res.print("res:");
#endif
    mat.ir(3,k).set(resp).rest(res).rs();
  }
  mat.scale(1/eps);
  x0.reshape(2,nel,ndim);
  mat.rs().reshape(4,nel,ndim,nel,ndim);

  if (0 && !glob_param->inwt) {
    x0.set(xloc).axpy(state_new,c_relax)
      .axpy(state_old,1.-c_relax);
    df_grad(x0,res);
    res.scale(1./c_relax);
  } 
  res.scale(-1.);

  dstate.set(state_new).rest(state_old);
  res_Dir.prod(mat,dstate,1,2,-1,-2,-1,-2);
  res.axpy(res_Dir,-1.);
  
#ifdef DEBUG_ANAL
  xloc.print("eig_anal: xloc");
  xloc.print("state_new");
  res.print("res:");
  mat.reshape(2,nel*ndim,nel*ndim).print("mat: ");
  printf("eps: %f\n",eps);
  PetscFinalize();
  exit(0);
#endif
}
