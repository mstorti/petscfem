/*
  This file belongs to he PETSc - FEM package a library and
  application suite oriented to the Finite Element Method based on PETSc. 
  Copyright (C) 1999, 2000  Mario Alberto Storti
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.

*/
  
#include "../../src/fem.h"
#include "../../src/utils.h"
#include "../../src/readmesh.h"
#include "../../src/getprop.h"

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

struct Swfm2tWallFunData {
  // Model of Launder & Spalding[1974], Krishnappan & Lau[1986]
  double von_Karman_cnst,roughness,A0,A1,A2,A3,A4,rho;
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "swfm2t_wall_res_fun" 
double swfm2t_wall_res_fun(double yn,double vwall,
			   double g,void *user_data) {

  Swfm2tWallFunData *s = (Swfm2tWallFunData *)user_data;
  double U_star = g/rho;
  double y_plus = U_star*yn/s->viscosity;
  double rough_nbr = U_star*s->roughness/s->viscosity;
  double logrou=log(rough_nbr);
  double exprou=exp(-s->A0*SQ(logrou));
  double Bs=exprou*(s->A1+s->A2*logrou)+s->A3*(1-exprou);
  double E=s->von_Karman_cnst*Bs/rough_nbr;
  double RHS=U_star/s->von_Karman_cnst*log(E*y_plus);
  return vwall-RHS;
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "swfm2t_wall_fun" 
double swfm2t_wall_res_fun(double yn,double vwall,
			   double &g, double &gprime,void *user_data) {
  // Secant iteration in order to find g(yn)
  Swfm2tWallFunData *s = (Swfm2tWallFunData *)user_data;
  double v1=vwall,g1=g;
  r1=swfm2t_wall_res_fun(yn,v1,g1,



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

  // Get arguments from arg_list

  //o The $y^+$ coordinate of the computational boundary
  NSGETOPTDEF(double,y_wall_plus,25.);
  double fwall,fprime;
  arg_data *staten,*stateo,*retval,*Ajac,*jac_prof;
  if (comp_res) {
    swfm2t_wall_fun(y_wall_plus,fwall,fprime);
    int j=-1;
    staten = &arg_data_v[++j];
    stateo = &arg_data_v[++j];
    retval  = &arg_data_v[++j];
    Ajac = &arg_data_v[++j];
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

  FastMat2 matloc(4,nel,ndof,nel,ndof);
  FastMat2 matlocmom(2,nel,nel);

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
  FMatrix grad_p_star(ndim),u,u_star,ucols,ucols_new,ucols_star,
    pcol_star,pcol_new,pcol;
  FMatrix matloc_prof(nen,nen),uc(ndim),tmp1,tmp2,tmp3,tmp4,tmp5,seed(ndof,ndof);

  if (comp_mat_res) {
    seed.eye().setel(0.,ndof,ndof);
  }

  if (comp_mat) {
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
      lambda_max=0;
      locstateo.set(element.vector_values(*stateo));
      locstaten.set(element.vector_values(*staten));
    }

    // State at time t_{n+\alpha}
    locstate.set(0.).axpy(locstaten,ALPHA).axpy(locstateo,(1-ALPHA));

    matlocmom.set(0.);
    matloc.set(0.);
    veccontr.set(0.);

    ucolso.set(locstateo.is(2,1,ndim));
    locstate2.rs();

    ucolsn.set(locstaten.is(2,1,ndim));
    locstate.rs();

    ucols_star.set(ucolsn).scale(alpha).axpy(ucolso,1-alpha);

    if (comp_res) {
      // loop over Gauss points
      veccontr.is(2,1,ndim);
      matloc.is(2,1,ndim).is(4,1,ndim);

      // Guarda que hay que invertir la direccion de la traccion!!!!
      for (ipg=0; ipg<npg; ipg++) {

	Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);
	detJaco = mydetsur(Jaco,normal);

	u_star.prod(SHAPE,ucols_star,-1,-1,1);
	double Ustar = sqrt(u_star.sum_square_all());
	double gfun,gprime;
	gprime = rho / (fwall*fwall);
	gfun = gprime * Ustar;

	tmp1.prod(SHAPE,SHAPE,1,2);
	matlocmom.axpy(tmp1,gfun*detJaco);

	tmp2.prod(SHAPE,u_star,1,2);
	tmp3.set(tmp2).scale(detJaco*gprime/Ustar);
	tmp4.prod(tmp2,tmp3,1,2,3,4);
	matloc.add(tmp4);

	veccontr.axpy(tmp2,-gfun*detJaco);
      }

      matlocmom.rs();
      tmp5.prod(matlocmom,seed,1,3,2,4);

      veccontr.rs();
      matloc.rs().add(tmp5);

      veccontr.export_vals(element.ret_vector_values(*retval));
#ifdef CHECK_JAC
      veccontr.export_vals(element.ret_fdj_values(*fdj_jac));
#endif
      veccontr.export_vals(&(RETVAL(ielh,0,0)));
      matloc.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
      matlocf.export_vals(element.ret_mat_values(*Ajac));
    } else if (comp_prof) {
      matlocf.export_vals(element.ret_mat_values(*jac_prof));
    }
    

  }
      
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  return 0;
}
#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
