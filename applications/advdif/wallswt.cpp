//__INSERT_LICENSE__
//$Id: wallswt.cpp,v 1.8 2001/05/25 04:12:29 mstorti Exp $
  
#include "../../src/fem.h"
#include "../../src/utils.h"
#include "../../src/readmesh.h"
#include "../../src/getprop.h"
#include "../../src/secant.h"

#include "../../src/fastmat2.h"
#include "advective.h"

extern TextHashTable *GLOBAL_OPTIONS;
extern int MY_RANK,SIZE;
extern int TSTEP; //debug:=
#define MAXPROP 100


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int wall::ask(char *,int &)"
int wall_swfm2t::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(comp_res);
  DONT_SKIP_JOBINFO(comp_prof);
  return 0;
}

class SwtWallFun : public Secant {
public:
  double y_wall_plus,von_Karman_cnst,roughness,A0,A1,A2,A3,rho,vwall,
    viscosity;
  SwtWallFun(double x0=0) : Secant(x0) {};
  double residual(double x,void * = NULL);
};


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "SwtWallFun::residual()"
double SwtWallFun::residual(double g,void *u=NULL) {

  double U_star = g/rho;
  // double y_plus = U_star*yn/viscosity;
  double rough_nbr = U_star * roughness / viscosity;
  double logrou=log(rough_nbr);
  double exprou=exp( - A0 * SQ(logrou) );
  double Bs=exprou*(A1 + A2 * logrou) + A3 * (1-exprou);
  double E=von_Karman_cnst * Bs / rough_nbr;
  double RHS=U_star / von_Karman_cnst * log(E * y_wall_plus);
  return vwall - RHS;
}



//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "wall::assemble"
void wall_swfm2t::new_assemble(arg_data_list &arg_data_v,const Nodedata *nodedata,
			      const Dofmap *dofmap,const char *jobinfo,
			      const ElementList &elemlist,
			      const TimeData *time_data) {

  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_prof);

  int ierr=0;

  NSGETOPTDEF(int,npg,0); //nd
  NSGETOPTDEF(int,ndim,0); //nd
  assert(npg>0);
  assert(ndim>0);

  int nel,ndof,nelprops;
  elem_params(nel,ndof,nelprops);
  int nen = nel*ndof;
  int ndimel = ndim-1;
  int nu=nodedata->nu;
  int nH = nu-ndim;

  // Get arguments from arg_list

  SwtWallFun swt_wall_fun;

#undef GETOPTDEF_HOOK
#define GETOPTDEF_HOOK(name) (swt_wall_fun.name)
  // Parameters for the wall law function
  //o Viscosity
  NSGETOPTDEF_ND(double,viscosity,1);
  //o The non-dmensional distance to the wall
  NSGETOPTDEF_ND(double,y_wall_plus,8.);
  //o The von Karman constant
  NSGETOPTDEF_ND(double,von_Karman_cnst,1.);
  //o Roughness of the wall
  NSGETOPTDEF_ND(double,roughness,0.);
  //o Parameter of the model
  NSGETOPTDEF_ND(double,A0,-0.217);
  //o Parameter of the model
  NSGETOPTDEF_ND(double,A1,5.50);
  //o Parameter of the model
  NSGETOPTDEF_ND(double,A2,2.50);
  //o Parameter of the model
  NSGETOPTDEF_ND(double,A3,8.5);
  // Return GETOPTDEF_ND_HOOK to his previous value
#undef GETOPTDEF_ND_HOOK
#define GETOPTDEF_ND_HOOK(name) name

  GlobParam *glob_param;
  arg_data *staten,*stateo,*retval,*Ajac,*jac_prof,*fdj_jac;
  if (comp_res) {
    int j=-1;
    staten = &arg_data_v[++j];
    stateo = &arg_data_v[++j];
    retval  = &arg_data_v[++j];
    Ajac = &arg_data_v[++j];
    ++j; // this is dtmin, we don't use this here
    glob_param = (GlobParam *) arg_data_v[++j].user_data;;
#ifdef CHECK_JAC
    fdj_jac = &arg_data_v[++j];
#endif
  }

  FastMat2 matlocf(4,nel,ndof,nel,ndof);
  if (comp_prof) {
    jac_prof = &arg_data_v[0];
    matlocf.set(1.);
  }

  // allocate local vecs
  int kdof;
  FMatrix veccontr(nel,ndof),xloc(nel,ndim),xc(ndim),locstate(nel,ndof),
    locstateo(nel,ndof),locstaten(nel,ndof); 

  FastMat2 matloc(4,nel,ndof,nel,ndof),Hloc(2,nel,nH);
  FastMat2 matlocmom(2,nel,nel);

  // The trapezoidal rule integration parameter 
#define ALPHA (glob_param->alpha)
#define DT (glob_param->Dt)

  double rho=1.;
  // Trapezoidal method parameter. 
  NSGETOPTDEF(double,alpha,1.);    //nd

  // Gauss Point data
  //o Type of element geometry to define Gauss Point data
  NGETOPTDEF_S(string,geometry,cartesian2d);
  
  GPdata gp_data(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);

  double detJaco,p_star,wpgdet;
  int elem, ipg,node, jdim, kloc,lloc,ldof;
  FMatrix Jaco(ndimel,ndim),resmom(nel,ndim),normal(ndim);
  FMatrix grad_p_star(ndim),u,u_star,ucolsn,ucolso,ucols_star,
    pcol_star,pcol_new,pcol;
  FMatrix matloc_prof(nen,nen),uc(ndim),tmp1,tmp2,tmp3,tmp4,tmp5,seed(ndof,ndof);

  if (comp_res) {
    seed.eye().setel(0.,ndof,ndof);
  }

  if (comp_prof) {
    matloc_prof.set(1.);
  }

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  for (ElementIterator element = elemlist.begin(); 
       element!=elemlist.end(); element++) {
    int k,ielh;
    element.position(k,ielh);

    FastMat2::reset_cache();

    element.node_data(nodedata,xloc.storage_begin(),
		      Hloc.storage_begin());

    if (comp_prof) {
      matlocf.export_vals(element.ret_mat_values(*jac_prof));
      continue;
    }

    if (comp_res) {
      locstateo.set(element.vector_values(*stateo));
      locstaten.set(element.vector_values(*staten));
    }

    // State at time t_{n+\alpha}
    locstate.set(0.).axpy(locstaten,ALPHA).axpy(locstateo,(1-ALPHA));

    matlocmom.set(0.);
    matloc.set(0.);
    veccontr.set(0.);

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    ucolso.set(locstateo.is(2,1,ndim));
    locstateo.rs();

    ucolsn.set(locstaten.is(2,1,ndim));
    locstaten.rs();

    ucols_star.set(ucolsn).scale(alpha).axpy(ucolso,1-alpha);

    if (comp_res) {
      // loop over Gauss points
      veccontr.is(2,1,ndim);
      matloc.is(2,1,ndim).is(4,1,ndim);

      // Guarda que hay que invertir la direccion de la traccion!!!!
      for (ipg=0; ipg<npg; ipg++) {

	Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);
	detJaco = mydetsur(Jaco,normal);
	normal.scale(-1.); // fixme:= This is to compensate a bug in mydetsur

	u_star.prod(SHAPE,ucols_star,-1,-1,1);
	double Ustar = sqrt(u_star.sum_square_all());
	swt_wall_fun.vwall = Ustar;
	double g = swt_wall_fun.sol();
#if 0
	gprime = rho / (fwall*fwall);
	gfun = gprime * Ustar;

	tmp1.prod(SHAPE,SHAPE,1,2);
	matlocmom.axpy(tmp1,gfun*detJaco);

	tmp2.prod(SHAPE,u_star,1,2);
	tmp3.set(tmp2).scale(detJaco*gprime/Ustar);
	tmp4.prod(tmp2,tmp3,1,2,3,4);
	matloc.add(tmp4);

	veccontr.axpy(tmp2,-gfun*detJaco);
#endif
      }

      matlocmom.rs();
      tmp5.prod(matlocmom,seed,1,3,2,4);

      veccontr.rs();
      matloc.rs().add(tmp5);

      veccontr.export_vals(element.ret_vector_values(*retval));
#ifdef CHECK_JAC
      veccontr.export_vals(element.ret_fdj_values(*fdj_jac));
#endif
      matlocf.export_vals(element.ret_mat_values(*Ajac));
    } else if (comp_prof) {
      matlocf.export_vals(element.ret_mat_values(*jac_prof));
    }
    

  }
      
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
}
#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
