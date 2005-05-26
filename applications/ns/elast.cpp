//__INSERT_LICENSE__
//$Id: elast.cpp,v 1.17 2005/05/26 22:07:29 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include "./fm2funm.h"

#include "nsi_tet.h"
#include "adaptor.h"
#include "elast.h"

class MyFun2 : public FastMat2_funm {
public:
  // Some kind of `hardening'
  // double f(double l) { return pow(l,-0.3); }

  // No hardening at all
  double f(double l) { return 1.; }
} my_fun2;

#if 0
class MyFun : public FastMat2_fund {
public:
  void f(const FastMat2 &L,FastMat2 &fL) { fL.set(L); }
} my_fun;
#endif

void elasticity::init() {

  int ierr;
  //o Young modulus
  TGETOPTDEF(thash,double,Young_modulus,0.);
  E=Young_modulus;
  assert(Young_modulus>0.);

  //o Poisson ratio
  TGETOPTDEF(thash,double,Poisson_ratio,0.);
  nu=Poisson_ratio;
  assert(nu>=0. && nu<0.5);

  //o Density
  TGETOPTDEF(thash,double,density,0.);
  rho=density;
  // Dos opciones para imprimir
  // printf("rec_Dt: %d\n",rec_Dt);
  // SHV(rec_Dt);
  assert(!(rec_Dt>0. && rho==0.));
  assert(ndof==ndim);

  //o exponent for stiffen deforming elements (Tezduyar paper)
  TGETOPTDEF(thash,double,stiffness_exponent,0.);
  Jaco_pow = stiffness_exponent;

  //o minimum det(Jaco) value for deformed mesh to stiffen deformed elements
  TGETOPTDEF(thash,double,detJaco_minimum_value,1.0e-10);
  detJaco_min = detJaco_minimum_value;

  G_body.resize(1,ndim).set(0.);
  const char *line;
  vector<double> G_body_v;
  thash->get_entry("G_body",line);
  if(line) {
    read_double_array(G_body_v,line);
    assert(G_body_v.size()==ndim);
    G_body.set(&G_body_v[0]);
  }

  ntens = ndim*(ndim+1)/2;
  nen=nel*ndof;
  
  // tal vez el resize blanquea
  B.resize(2,ntens,nen).set(0.);
  C.resize(2,ntens,ntens).set(0.);
  Jaco.resize(2,ndim,ndim);
  Jaco_def.resize(2,ndim,ndim);
  G.resize(2,ndim,ndim);
  dshapex.resize(2,ndim,nel);  
  dJaco.resize(2,nel,ndim);
  tmp_elast.resize(1,nel);

  my_fun2.init(G);
  // Plane strain
  if (ndim==2) {
    double c1=E*(1.-nu)/((1.+nu)*(1.-2.*nu)), c2=E/(2.*(1.+nu)),
      c3=nu/(1.-nu);
    C.setel(c1,1,1)
      .setel(c1*c3,1,2)
      .setel(c1*c3,2,1)
      .setel(c1,2,2)
      .setel(c2,3,3);
  } else if (ndim==3) {
    double c1=E*(1.-nu)/((1.+nu)*(1.-2.*nu)), 
      c2 = (1-2*nu)/2./(1-nu),
      c3=nu/(1.-nu);
      C.is(1,1,3).is(2,1,3).set(c3)
	.rs().d(1,2).is(1,1,3).set(1.)
	.rs().d(1,2).is(1,4,6).set(c2)
	.rs().scale(c1);
  } else {
    PetscPrintf(PETSC_COMM_WORLD,"wrong dimension: %d\n",ndim);
    assert(0);
  }

}

void elasticity::element_connector(const FastMat2 &xloc,
				   const FastMat2 &state_old,
				   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){

  B.reshape(3,ntens,nel,ndof);
  
  // Levi-Civita tensor generation
  epsilon_LC.eps_LC();
  double J1,J2,coef;
  
  // loop over Gauss points
  for (int ipg=0; ipg<npg; ipg++) {
    
    x_new.set(xloc);
    //**    x_new.set(xloc).add(state_new);
    dshapexi.ir(3,ipg+1); // restriccion del indice 3 a ipg+1
    //x_def.set(xloc).add(state_old);
    x_def.set(xloc).add(state_new);
    Jaco_def.prod(dshapexi,x_def,1,-1,-1,2);
    Jaco.prod(dshapexi,x_new,1,-1,-1,2);
    G.prod(Jaco,Jaco,-1,1,-1,2);
    my_fun2.apply(G,fG);
    
    double detJaco = Jaco.det();
    double detJaco_def = Jaco_def.det();
    //**    if (detJaco_0<=0.) {
    if (detJaco<=0.) {
      detj_error(detJaco,elem);
      set_error(1);
    }
    detJaco_def = (detJaco_def>detJaco_min ? detJaco_def : detJaco_min);
    double wpgdet = detJaco*wpg.get(ipg+1);
    wpgdet = wpgdet*pow(detJaco/detJaco_def,Jaco_pow);
    iJaco.inv(Jaco);
    dshapex.prod(iJaco,dshapexi,1,-1,-1,2);
    dshapex_scaled.prod(fG,dshapex,1,-1,-1,2);
    //dshapex_scaled.set(dshapex);
    
    // Recall: \epsilon = B dudx
    // where \epsilon = [e_xx e_yy e_xy] (2D)
    //                  [e_xx e_yy e_zz e_xy e_xz e_yz] (3D)
    // dudx is de gradient of displacements. 
    if (ndim==2) {
      B.ir(1,1).ir(3,1).set(dshapex_scaled.ir(1,1));
      B.ir(1,2).ir(3,2).set(dshapex_scaled.ir(1,2));
      B.ir(1,3).ir(3,1).set(dshapex_scaled.ir(1,2));
      B.ir(1,3).ir(3,2).set(dshapex_scaled.ir(1,1));
    } else if (ndim==3) {
      B.ir(1,1).ir(3,1).set(dshapex_scaled.ir(1,1));
      B.ir(1,2).ir(3,2).set(dshapex_scaled.ir(1,2));
      B.ir(1,3).ir(3,3).set(dshapex_scaled.ir(1,3));

      B.ir(1,4).ir(3,2).set(dshapex_scaled.ir(1,1));
      B.ir(1,4).ir(3,1).set(dshapex_scaled.ir(1,2));

      B.ir(1,5).ir(3,3).set(dshapex_scaled.ir(1,1));
      B.ir(1,5).ir(3,1).set(dshapex_scaled.ir(1,3));

      B.ir(1,6).ir(3,3).set(dshapex_scaled.ir(1,2));
      B.ir(1,6).ir(3,2).set(dshapex_scaled.ir(1,3));
    } else assert(0);

    dshapex_scaled.rs();
    
    // B.rs().reshape(2,ntens,nen);
    B.rs();
    strain.prod(B,state_new,1,-1,-2,-1,-2);
    stress.prod(C,strain,1,-1,-1);
    
    // Residual computation
    res_pg.prod(B,stress,-1,1,2,-1);
    res.axpy(res_pg,-wpgdet);

    // Jacobian computation
    mat_pg1.prod(C,B,1,-1,-1,2,3);
    mat_pg2.prod(B,mat_pg1,-1,1,2,-1,3,4);
    mat.axpy(mat_pg2,wpgdet);
    
    dJaco.set(0.);
    if (ndim==2) {
      J1 = Jaco_def.get(2,2);
      dshapexi.ir(1,1);
      tmp_elast.set(dshapexi).scale(J1);
      dJaco.ir(2,1).set(tmp_elast).rs();

      J1 = Jaco_def.get(1,2);
      dshapexi.ir(1,2);
      tmp_elast.set(dshapexi).scale(-J1);
      dJaco.ir(2,1).add(tmp_elast).rs();

      J1 = Jaco_def.get(1,1);
      dshapexi.ir(1,2);
      tmp_elast.set(dshapexi).scale(J1);
      dJaco.ir(2,2).set(tmp_elast).rs();

      J1 = Jaco_def.get(2,1);
      dshapexi.ir(1,1);
      tmp_elast.set(dshapexi).scale(-J1);
      dJaco.ir(2,2).add(tmp_elast).rs();
      
    } else if (ndim==3) {
      for (int p=1; p<=3; p++) {
	for (int q=1; q<=3; q++) {
	  for (int r=1; r<=3; r++) {
	    //epsilon_LC.ir(1,p).ir(2,q).ir(3,r);
	    coef = epsilon_LC.get(p,q,r);
	    // first term
	    J1 = Jaco_def.get(q,2);
	    J2 = Jaco_def.get(r,3);
	    dshapexi.ir(1,p);
	    //tmp_elast.set(dshapexi).scale(J1*J2*double(epsilon_LC));
	    tmp_elast.set(dshapexi).scale(J1*J2*coef);
            dJaco.ir(2,1).add(tmp_elast).rs();

	    // second term
	    J1 = Jaco_def.get(p,1);
	    J2 = Jaco_def.get(r,3);
	    dshapexi.ir(1,q);
	    //	    tmp_elast.set(dshapexi).scale(J1*J2*double(epsilon_LC));
	    tmp_elast.set(dshapexi).scale(J1*J2*coef);
	    dJaco.ir(2,2).add(tmp_elast).rs();

	    // third term
	    J1 = Jaco_def.get(p,1);
	    J2 = Jaco_def.get(q,2);
	    dshapexi.ir(1,r);
	    //	    tmp_elast.set(dshapexi).scale(J1*J2*double(epsilon_LC));
	    tmp_elast.set(dshapexi).scale(J1*J2*coef);
	    dJaco.ir(2,3).add(tmp_elast).rs();
	    
	    //	    epsilon_LC.rs();

	  }
	}
      }
      //      dshapexi.rs();
    }
    dshapexi.rs();

    mat_pg3.prod(res_pg,dJaco,1,2,3,4).scale(-Jaco_pow/detJaco_def);
    mat.axpy(mat_pg3,wpgdet);
  }
    
}

