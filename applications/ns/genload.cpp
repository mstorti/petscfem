//__INSERT_LICENSE__
//$Id: genload.cpp,v 1.6 2002/12/16 22:28:24 mstorti Exp $
#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/util2.h>
#include <src/fastmat2.h>

#include "./nsi_tet.h"
#include "./genload.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int ns_volume_element::ask(char *,int &)"
int GenLoad::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(comp_mat);
  DONT_SKIP_JOBINFO(comp_res);
  DONT_SKIP_JOBINFO(comp_mat_res);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void GenLoad::q(FastMat2 &uin,FastMat2 &uout,FastMat2 &flux,
		FastMat2 &jacin,FastMat2 &jacout) {
  PETSCFEM_ERROR0("Not defined one layer flux function!\n");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void GenLoad::q(FastMat2 &uin,FastMat2 &flux,FastMat2 &jacin) {
  PETSCFEM_ERROR0("Not defined double layer flux function!\n");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "GenLoad::assemble"
int GenLoad::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
		      Dofmap *dofmap,const char *jobinfo,int myrank,
		      int el_start,int el_last,int iter_mode,
		      const TimeData *time_) {

  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_mat_res);
  GET_JOBINFO_FLAG(comp_res);

  // Essentially treat comp_res as comp_mat_res but
  // with the side effect of update_jacobian=1
  int update_jacobian=1;	
  if (comp_res) {
    comp_mat_res=1;
    update_jacobian=0;
  }

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retvalmat,iele,j,nel,k,ndof,p,nel,q,ndof)
#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 

  int ierr;
  TGETOPTNDEF(thash,int,ndim,none); //nd
  int nen = nel*ndof;

  // Get arguments from arg_list
  double *locst,*locst2,*retval,*retvalmat;
  if (comp_mat) {
    retvalmat = arg_data_v[0].retval;
  } 

  // rec_Dt is the reciprocal of Dt (i.e. 1/Dt)
  // for steady solutions it is set to 0. (Dt=inf)
  GlobParam *glob_param;
  double *hmin,Dt,rec_Dt;
  int ja_hmin;
#define WAS_SET arg_data_v[ja_hmin].was_set
  if (comp_mat_res) {
    int ja=0;
    locst = arg_data_v[ja++].locst;
    locst2 = arg_data_v[ja++].locst;
    retval = arg_data_v[ja++].retval;
    if (update_jacobian) retvalmat = arg_data_v[ja++].retval;
    hmin = (arg_data_v[ja++].vector_assoc)->begin();
    ja_hmin=ja;
    glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
    rec_Dt = 1./glob_param->Dt;
    if (glob_param->steady) rec_Dt=0.;
  } 
  double &alpha = glob_param->alpha;

  int ndimel = ndim-1;
  // Unpack nodedata
  int nu=nodedata->nu;
  //o Number of Gauss points.
  TGETOPTNDEF(thash,int,npg,none);
  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  //o Whether there is a double or single layer of nodes
  TGETOPTDEF(thash,int,double_layer,0);

  nel2;
  if (double_layer) {
    PETSCFEM_ASSERT0(nel % 2 ==0,"Number of nodes per element has to be even for "
		    "double_layer mode");
    nel2=nel/2;
  } else {
    nel2=nel;
  }

  assert(nel % 2==0); // one layer of nodes is for \eta the other for w
  int nH = 0;

  GPdata gp_data(geometry.c_str(),ndimel,nel2,npg,GP_FASTMAT2);

  FastMat2 matloc_prof(4,nel,ndof,nel,ndof),
    matlocf(4,nel,ndof,nel,ndof),
    state_star(2,nel,ndof), state_old(2,nel,ndof),
    u_in(2,nel2,ndof), u_out(2,nel2,ndof), 
    U_in(1,ndof), U_out(1,ndof), 
    flux_in(1,ndof), flux_out(1,ndof), 
    jac(2,ndof,ndof), veccontr(2,nel,ndof), 
    Jaco, iJaco, dshapex(2,ndimel,nel2), xloc(2,nel2,ndim),
    h_in, h_out, flux, jac_in, jac_out, tmp1,
    tmp2, tmp3, tmp4, vecc2;

  if (double_layer) {
    jac.resize(4,2,2,ndof,ndof);
  }

  veccontr.set(0.);
  matlocf.set(0.);

  if (comp_mat) matloc_prof.set(1.);
  start_chunk();

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  int ielh=-1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    FastMat2::reset_cache();
    ielh++;
    // load_props(propel,elprpsindx,nprops,&(ELEMPROPS(k,0)));
    int elem = k;
    element_hook(k);

    // Load local node coordinates in local vector
    for (int kloc=0; kloc<nel2; kloc++) {
      int node = ICONE(k,kloc);
      xloc.ir(1,kloc+1).set(&NODEDATA(node-1,0));
    }
    xloc.rs();

    matlocf.set(0.);
    veccontr.set(0.);

    if(comp_mat) {
      matloc_prof.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
      continue;
    }

    state_old.set(&(LOCST2(ielh,0,0)));
    state_star.set(&(LOCST(ielh,0,0))).scale(alpha).axpy(state_old,1.-alpha);
    state_star.is(1,1,nel2);
    u_in.set(state_star);
    if (double_layer) {
      state_star.rs().is(1,nel2+1,nel);
      u_out.set(state_star);
    }

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    for (int ipg=0; ipg<npg; ipg++) {

      Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);
      double detJaco = Jaco.detsur();
      if (detJaco <= 0.) {
	printf("Jacobian of element %d is negative or null\n"
	       " Jacobian: %f\n",k,detJaco);
	PetscFinalize();
	exit(0);
      }
      double wpgdet = detJaco*WPG;

      U_in.prod(SHAPE,u_in,-1,-1,1);
      if (nH>0) H_m.prod(SHAPE,h_in,-1,-1,1);
      if (double_layer) {
	U_out.prod(SHAPE,u_out,-1,-1,1);
	if (nH>0) H_out_m.prod(SHAPE,h_out,-1,-1,1);
	q(U_in,U_out,flux_in,flux_out,jac);
      } else {
	assert(0);
	// q(U_in,flux,jac_in);
      }
      
      tmp1.set(SHAPE).scale(wpgdet);
      vecc2.prod(tmp1,flux_in,1,2);

      veccontr.is(1,1,nel2).add(vecc2);
      veccontr.rs();
      // Contribution to jacobian from interior side
      tmp2.set(SHAPE).scale(wpgdet);
      tmp3.prod(SHAPE,tmp2,1,2);
      jac.ir(1,1).ir(2,1);
      tmp4.prod(tmp3,jac,1,3,2,4);
      matlocf.is(1,1,nel2).is(3,1,nel2).add(tmp4);
      jac.rs();
      if (double_layer) {
	matlocf.reshape(6,2,nel2,ndof,2,nel2,ndof);
	for (int b=1; b<=2; b++) {
	  for (int bb=1; bb<=2; bb++) {
	    if (b==1 && bb==1) continue;
	    tmp2.set(SHAPE).scale(wpgdet);
	    tmp3.prod(SHAPE,tmp2,1,2);
	    jac.ir(1,b).ir(2,bb);
	    tmp4.prod(tmp3,jac,1,3,2,4);
	    matlocf.ir(1,b).ir(4,bb);
	    matlocf.add(tmp4);
	  }
	}
	matlocf.rs(); jac.rs();
	vecc2.prod(tmp1,flux_out,1,2);
	veccontr.is(1,nel2+1,nel).add(vecc2);
	veccontr.rs();
	matlocf.reshape(4,nel,ndof,nel,ndof);
      }
      matlocf.rs();
    }
    veccontr.export_vals(&(RETVAL(ielh,0,0)));
    matlocf.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  return 0;
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void ConsGenLoad::q(FastMat2 &uin,FastMat2 &uout,
		    FastMat2 &flux_in, FastMat2 &flux_out, 
		    FastMat2 &jac) {
  q(uin,uout,flux_in,jac_aux);
  flux_out.set(flux_in).scale(-1.);

  // Jac(2,:) = -Jac(1,:)
  jac.ir(1,1).set(jac_aux)
    .ir(1,2).set(jac_aux).scale(-1).rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void ConsGenLoad::start_chunk() {
  jac_aux.resize(3,2,ndof,ndof);
  start_chunk_c();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void lin_gen_load::start_chunk_c() {
  int ierr;
  //o The film coefficient constant
  TGETOPTDEF_ND(thash,double,h_film,0);
  assert(h_film>0.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void lin_gen_load::q(FastMat2 &U_in,FastMat2 &U_out,
		     FastMat2 &flux_in,FastMat2 &jac) {
  flux_in.set(U_out).rest(U_in).scale(h_film);
  jac.ir(1,1).eye(h_film).rs().ir(1,2).eye(-h_film).rs();
}

