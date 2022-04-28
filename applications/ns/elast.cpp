//__INSERT_LICENSE__
//$Id: elast.cpp,v 1.22.10.1 2007/02/19 20:23:56 mstorti Exp $

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
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define MAXPROPS 100
  elprpsindx.mono(MAXPROPS);
  propel.mono(MAXPROPS);
  
  // TGETOPTNDEF_ND(thash,int,ndim,none); //nd

  int iprop=0;
  Young_modulus_indx = iprop; 
  ierr = get_prop(iprop,elem_prop_names,
		  thash,elprpsindx.buff(),propel.buff(), 
		  "Young_modulus",1);
  nprops = iprop;

  //o Scale Young modulus
  TGETOPTDEF_ND(thash,double,Young_modulus_fac,1.0);
  assert(Young_modulus_fac>0.);

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
    assert(G_body_v.size()==(unsigned int)ndim);
    G_body.set(&G_body_v[0]);
  }

  ntens = ndim*(ndim+1)/2;
  nen=nel*ndof;

  // Dump elemset values
  TGETOPTDEF_ND(thash,int,per_element_vals_size,0);
  int dump_elem_vals = (per_element_vals_size>0);
  if (dump_elem_vals) {
    PETSCFEM_ASSERT0(per_element_vals_size==ntens,
                     "elasticity element stores stress values only");
    // FIXME:= this fails
    if (!!per_element_vals_p) {
      per_element_vals_p.reset(new dvector<double>);
      dvector<double> &vals = *per_element_vals_p;
      printf("vals %p\n",&vals);
      vals.a_resize(2,nelem,ntens);
      vals.defrag();
      vals.set(0.0);
      exit(0);
    }
  }
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
    double c1=(1.-nu)/((1.+nu)*(1.-2.*nu)), c2=1.0/(2.*(1.+nu)),
      c3=nu/(1.-nu);
    C.setel(c1,1,1)
      .setel(c1*c3,1,2)
      .setel(c1*c3,2,1)
      .setel(c1,2,2)
      .setel(c2,3,3);
#if 0
    // From RAPPORT SF-101, LTAS, Godinas A., Jetteur P., Laschet G.
    // This modification of the constitutive relation should avoid
    // blockage of thin elements for plate model. 
    int plate_model = 1;
    if (plate_model) {
      printf("using plate model\n");
      C.set(0.0);
      double c1 = 1.0/(1.0-nu*nu);
      C.setel(c1,1,1);
      C.setel(1.0,2,2);
      C.setel(c1*(5.0/6.0)*(1.0-nu)/2.0,3,3);
    }
#endif
  } else if (ndim==3) {
    double c1=(1.-nu)/((1.+nu)*(1.-2.*nu)), 
      c2 = (1-2*nu)/2./(1-nu),
      c3=nu/(1.-nu);
      C.is(1,1,3).is(2,1,3).set(c3)
	.rs().d(1,2).is(1,1,3).set(1.)
	.rs().d(1,2).is(1,4,6).set(c2)
	.rs().scale(c1);
  } else {
    PetscPrintf(PETSCFEM_COMM_WORLD,"wrong dimension: %d\n",ndim);
    assert(0);
  }

}

void elasticity::element_connector(const FastMat2 &xloc,
				   const FastMat2 &state_old,
				   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){

  load_props(propel.buff(),elprpsindx.buff(),nprops,
	     &(ELEMPROPS(elem,0)));
  double Young_modulus = 
    *(propel.buff()+Young_modulus_indx)
    *Young_modulus_fac;
#if 0
  if (rand()%100==0) 
    printf("Young_modulus: %f, fac %f\n",
	   Young_modulus,Young_modulus_fac);
#endif

  B.reshape(3,ntens,nel,ndof);
  
  // Levi-Civita tensor generation
  epsilon_LC.eps_LC();
  double J1,J2,coef;
  int dump_elem_vals = (per_element_vals_size>0);
  
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
    stress.prod(C,strain,1,-1,-1).scale(Young_modulus);
    if (elem%50==0) {
      SHV(elem);
      FMSHV(stress);
    }
    if (dump_elem_vals) {
      dvector<double> &vals = *per_element_vals_p;
      // Store the stress for postprocessing
      double *w = stress.storage_begin();
      for (int j=0; j<ntens; j++) {
        PETSCFEM_ASSERT0((elem*ntens+j)<int(vals.size()),
                         "bad size");
        vals.e(elem,j) += w[j];
      }
    }
    
    // Residual computation
    res_pg.prod(B,stress,-1,1,2,-1);
    res.axpy(res_pg,-wpgdet);

    // Jacobian computation
    mat_pg1.prod(C,B,1,-1,-1,2,3).scale(Young_modulus);
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
  if (dump_elem_vals) {
    // Store the stress for postprocessing
    dvector<double> &vals = *per_element_vals_p;
    for (int j=0; j<ntens; j++) vals.e(elem,j) /= ntens;
    if (elem%50==0) {
      SHV(elem);
      for (int j=0; j<ntens; j++)
        printf("%f ",vals.e(elem,j));
      printf("\n");
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void elasticity::clean() {
  printf("in elasticity::clean()\n");
}
