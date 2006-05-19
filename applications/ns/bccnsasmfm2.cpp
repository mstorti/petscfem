//__INSERT_LICENSE__
//$Id: bccnsasmfm2.cpp,v 1.2.58.1 2006/05/19 23:43:37 dalcinl Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>

#include <src/fastmat2.h>
#include "nsi_tet.h"

extern TextHashTable *GLOBAL_OPTIONS;
#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "int bcconv_nsasm_fm2::ask(char *,int &)"
int bcconv_nsasm_fm2::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(comp_mat);
  DONT_SKIP_JOBINFO(comp_res);
  DONT_SKIP_JOBINFO(comp_mat_res);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "bcconv_nsasm_fm2::assemble"
int bcconv_nsasm_fm2::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
			  Dofmap *dofmap,const char *jobinfo,int myrank,
			  int el_start,int el_last,int iter_mode,
			  const TimeData *time_data) {

  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_mat_res);
  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(build_nneighbor_tree);

// added for scalar transport equation
  GET_JOBINFO_FLAG(comp_mat_th);
  GET_JOBINFO_FLAG(comp_mat_res_th);
  GET_JOBINFO_FLAG(comp_res_th);

  comp_mat_res_th=comp_mat_res;
  comp_mat_th=comp_mat;

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retvalmat,iele,j,nel,k,ndof,p,nel,q,ndof)

  int ierr=0;

#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)])
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define IDENT(j,k) (ident[ndof*(j)+(k)])
#define JDOFLOC(j,k) VEC2(jdofloc,j,k,ndof)

  int locdof,kldof,lldof;
  char *value;

  // Unpack Elemset
  int npg,ndim;
  ierr = get_int(thash,"npg",&npg); CHKERRA(ierr);
  ierr = get_int(thash,"ndim",&ndim); CHKERRA(ierr);

  int ndimel = ndim-1;
  int nen = nel*ndof;

  // Unpack Dofmap
  int *ident,neq,nnod;
  neq = dofmap->neq;
  nnod = dofmap->nnod;

  // Unpack nodedata
  int nu=nodedata->nu;
  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  // Get arguments from arg_list
  double *locst,*locst2,*retval,*retvalmat;
  if (comp_mat) {
    retvalmat = arg_data_v[0].retval;
  }

  if (comp_mat_res || comp_mat_res_th) {
    locst = arg_data_v[0].locst;
    locst2 = arg_data_v[1].locst;
    retval = arg_data_v[2].retval;
    retvalmat = arg_data_v[3].retval;
  }

  //o Use a weak form for the gradient of pressure term and transport equation
  SGETOPTDEF(int,weak_form,1);

  //o Density
  SGETOPTDEF(double,rho,1.);
  //o Specific heat - constant pressure
  SGETOPTDEF(double,Cp,1.);
  //% Bubble radius
  SGETOPTDEF(double,rb,0.001);
  //% gravity direction
  SGETOPTDEF(int,g_dir,3);

  Cp = Cp/rho;

  // allocate local vecs
  int kdof;
  FMatrix veccontr(nel,ndof),xloc(nel,ndim),locstate(nel,ndof),
    locstate2(nel,ndof),xpg;

  if (ndof != ndim+2) {
    PetscPrintf(PETSCFEM_COMM_WORLD,"ndof != ndim+2\n"); CHKERRA(1);
  }

  nen = nel*ndof;
  FMatrix matloc(nen,nen), matlocmom(nel,nel), matlocther(nel,nel);

  // Physical properties
  int iprop=0, elprpsindx[MAXPROP]; double propel[MAXPROP];

  // Trapezoidal method parameter.
  TGETOPTDEF(GLOBAL_OPTIONS,double,alpha,1.); //nd

  double pi = 4*atan(1.0);

  int nprops=iprop;

  // double rho=1.;

  // Gauss Point data
  //    char *geom;
  //    thash->get_entry("geometry",geom);

  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  GPdata gp_data(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);

  // Definiciones para descargar el lazo interno
  double detJaco,p_star,T_star,wpgdet;

  int elem, ipg,node, jdim, kloc,lloc,ldof;

  FMatrix Jaco(ndimel,ndim),resmom(nel,ndim),normal(ndim),matij(ndof,ndof),resther(nel);

  FMatrix grad_p_star(ndim),u,u_star,ucols,ucols_new,ucols_star,
    pcol_star,pcol_new,pcol,tmp1,tmp2,tmp3,tmp4,tmp2_th;

// modif A-C
  FMatrix u_gas_star;
// ==========================

  FMatrix Tcol_star,Tcol_new,Tcol;

  FMatrix matloc_prof(nen,nen);

  double tmp3_th,tmp4_th;

  if (comp_mat || comp_mat_th)  {
    matloc_prof.set(1.);
  }

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  int ielh=-1;
  int start_chunk=1; // BETO : esto parece que no se usa

  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    FastMat2::reset_cache();
    ielh++;
    load_props(propel,elprpsindx,nprops,&(ELEMPROPS(k,0)));

    // Load local node coordinates in local vector
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(k,kloc);
      xloc.ir(1,kloc+1).set(&NODEDATA(node-1,0));
    }
    xloc.rs();

    // tenemos el estado locstate2 <- u^n
    //                   locstate  <- u^*
    if (comp_mat_res || comp_mat_res_th) {
      locstate.set(&(LOCST(ielh,0,0)));
      locstate2.set(&(LOCST2(ielh,0,0)));
    }

    matlocmom.set(0.);
    matlocther.set(0.);
    matloc.set(0.);
    veccontr.set(0.);
    resmom.set(0.);
    resther.set(0.);

    ucols.set(locstate2.is(2,1,ndim));
    pcol.set(locstate2.rs().ir(2,ndim+1));
    Tcol.set(locstate2.rs().ir(2,ndof));
    locstate2.rs();

    ucols_new.set(locstate.is(2,1,ndim));
    pcol_new.set(locstate.rs().ir(2,ndim+1));
    Tcol_new.set(locstate.rs().ir(2,ndof));
    locstate.rs();

    ucols_star.set(ucols_new).scale(alpha).axpy(ucols,1-alpha);
    pcol_star.set(pcol_new).scale(alpha).axpy(pcol,1-alpha);
    Tcol_star.set(Tcol_new).scale(alpha).axpy(Tcol,1-alpha);

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    matij.set(0.);

    // loop over Gauss points
    for (ipg=0; ipg<npg; ipg++) {

      Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);
      detJaco = Jaco.detsur(&normal);
      if (detJaco<=0.) {
	detj_error(detJaco,elem);
	set_error(1);
      }

      // WPG is used instead of wpgdet because normal is not a unit vector ;

      if (comp_mat_res) {

	// state variables and gradient
	p_star = double(tmp1.prod(SHAPE,pcol_star,-1,-1));
	u_star.prod(SHAPE,ucols_star,-1,-1,1);

	// Residual contribution due to integration by parts of pressure gradient
	tmp2.set(normal).scale(p_star);
	tmp3.prod(SHAPE,tmp2,1,2);
	resmom.axpy(tmp3,-WPG);

        matij.set(0.);
	// Matrix contribution due to integration by parts of pressure gradient
	for (int iloc=1; iloc<=nel; iloc++) {
	  double SHAPE_iloc = SHAPE.get(iloc);
	  for (int jloc=1; jloc<=nel; jloc++) {
	    double SHAPE_jloc = SHAPE.get(jloc);
	    tmp2.set(normal).scale(SHAPE_jloc);
	    tmp2.scale(SHAPE_iloc);
	    matij.ir(2,ndim+1).is(1,1,ndim).set(tmp2);
	    matij.rs();

	    int il1=(iloc-1)*ndof+1;
	    int il2=il1+ndof-1;
	    int jl1=(jloc-1)*ndof+1;
	    int jl2=jl1+ndof-1;
	    matloc.is(1,il1,il2).is(2,jl1,jl2).axpy(matij,WPG).rs();

	  }
	}

      } else if (comp_mat) {
	// don't make anything here !!
      } else {
	PetscPrintf(PETSCFEM_COMM_WORLD,
		    "Don't know how to compute jobinfo: %s\n",jobinfo);
	CHKERRQ(ierr);
      }

      if (comp_mat_res_th) {
	// state variables and gradient

// modif
    double vslip = (rb<7e-4 ? 4474*pow(rb,1.357) :
  		  rb<5.1e-3 ? 0.23 : 4.202*pow(rb,0.547));
// ==========================

	u_star.prod(SHAPE,ucols_star,-1,-1,1);
	T_star = double(tmp1.prod(SHAPE,Tcol_star,-1,-1));

// modif
    u_gas_star.set(u_star).addel(vslip,g_dir);
// ==========================

	// Residual contribution due to integration by parts of thermal convection term
// modif
        tmp3_th = rho*Cp*double(tmp2_th.prod(normal,u_gas_star,-1,-1));
// ==========================
	tmp4_th = T_star*tmp3_th;
	resther.axpy(SHAPE,tmp4_th*(-WPG));

	// Matrix contribution due to integration by parts of thermal convective flux
        matij.set(0.);
	for (int iloc=1; iloc<=nel; iloc++) {
	  double SHAPE_iloc = SHAPE.get(iloc);
	  for (int jloc=1; jloc<=nel; jloc++) {
	    double SHAPE_jloc = SHAPE.get(jloc);
	    tmp4_th = SHAPE_iloc*SHAPE_jloc*tmp3_th;

	    matij.ir(2,ndof).ir(1,ndof).set(tmp4_th);
	    matij.rs();

	    int il1=(iloc-1)*ndof+1;
	    int il2=il1+ndof-1;
	    int jl1=(jloc-1)*ndof+1;
	    int jl2=jl1+ndof-1;
	    matloc.is(1,il1,il2).is(2,jl1,jl2).axpy(matij,WPG).rs();

	  }
	}

      } else if (comp_mat_th) {
	// don't make anything here !!
      } else {
	PetscPrintf(PETSCFEM_COMM_WORLD,
		    "Don't know how to compute jobinfo: %s\n",jobinfo);
	CHKERRQ(ierr);
      }

    }

    // BETO : matloc_prof habria que inicializarlo a cero ???????
    if(comp_mat) {
      matloc_prof.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
    }

    if (comp_mat_res || comp_mat_res_th) {
      if (!weak_form) {
	veccontr.set(0.);
	matloc.set(0.);
      } else {
        veccontr.is(2,1,ndim).set(resmom).rs();
        veccontr.ir(2,ndof).set(resther).rs();
      }

      veccontr.export_vals(&(RETVAL(ielh,0,0)));
      matloc.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
    }
  }

  FastMat2::void_cache();
  FastMat2::deactivate_cache();

}

#undef SHAPE
#undef DSHAPEXI
#undef WPG
