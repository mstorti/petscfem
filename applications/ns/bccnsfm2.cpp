//__INSERT_LICENSE__
//$Id: bccnsfm2.cpp,v 1.9 2001/12/20 21:58:55 mstorti Exp $
  
#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>

#include <src/fastmat2.h>
#include "nsi_tet.h"

// BETO : esto se usa ???  #include <src/util2.h>
// BETO : esto se usa ???  #define USE_FASTMAT

// BETO : esto se usa ???  
extern TextHashTable *GLOBAL_OPTIONS;
#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int bcconv_ns_fm2::ask(char *,int &)"
int bcconv_ns_fm2::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(comp_mat);
  DONT_SKIP_JOBINFO(comp_res);
  DONT_SKIP_JOBINFO(comp_mat_res);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "bcconv_ns_fm2::assemble"
int bcconv_ns_fm2::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
			  Dofmap *dofmap,const char *jobinfo,int myrank,
			  int el_start,int el_last,int iter_mode,
			  const TimeData *time_data) {

  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_mat_res);
  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(build_nneighbor_tree);

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

  if (comp_mat_res) {
    locst = arg_data_v[0].locst;
    locst2 = arg_data_v[1].locst;
    retval = arg_data_v[2].retval;
    retvalmat = arg_data_v[3].retval;
  }

  //o Use a weak form for the gradient of pressure term. 
  SGETOPTDEF(int,weak_form,1);

  // allocate local vecs
  int kdof;
  FMatrix veccontr(nel,ndof),xloc(nel,ndim),locstate(nel,ndof),
    locstate2(nel,ndof),xpg; 

  if (ndof != ndim+1) {
    PetscPrintf(PETSC_COMM_WORLD,"ndof != ndim+1\n"); CHKERRA(1);
  }

  nen = nel*ndof;
  FMatrix matloc(nen,nen), matlocmom(nel,nel);

  // Physical properties
  int iprop=0, elprpsindx[MAXPROP]; double propel[MAXPROP];

  // Trapezoidal method parameter. 
  TGETOPTDEF(GLOBAL_OPTIONS,double,alpha,1.); //nd

  double pi = 4*atan(1.0);

  int nprops=iprop;
  
  double rho=1.;

  // Gauss Point data
  //    char *geom;
  //    thash->get_entry("geometry",geom);
  
  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  GPdata gp_data(geometry.c_str(),ndim,nel,npg,GP_FASTMAT2);

  // Definiciones para descargar el lazo interno
  double detJaco,p_star,wpgdet;

  int elem, ipg,node, jdim, kloc,lloc,ldof;
    
  FMatrix Jaco(ndimel,ndim),resmom(nel,ndim),normal(ndim),matij(ndof,ndof);
          
  FMatrix grad_p_star(ndim),u,u_star,ucols,ucols_new,ucols_star,
    pcol_star,pcol_new,pcol,tmp1,tmp2,tmp3,tmp4;

  FMatrix matloc_prof(nen,nen);

  if (comp_mat) {
    matloc_prof.set(1.);
  }
    
#ifdef USE_FASTMAT2_CACHE
  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);
#endif

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
    if (comp_mat_res) {
      locstate.set(&(LOCST(ielh,0,0)));
      locstate2.set(&(LOCST2(ielh,0,0)));
    }

    matlocmom.set(0.);
    matloc.set(0.);
    veccontr.set(0.);
    resmom.set(0.);

    ucols.set(locstate2.is(2,1,ndim));
    pcol.set(locstate2.rs().ir(2,ndof));
    locstate2.rs();

    ucols_new.set(locstate.is(2,1,ndim));
    pcol_new.set(locstate.rs().ir(2,ndof));
    locstate.rs();

    ucols_star.set(ucols_new).scale(alpha).axpy(ucols,1-alpha);
    pcol_star.set(pcol_new).scale(alpha).axpy(pcol,1-alpha);

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    matij.set(0.);

    // loop over Gauss points
    for (ipg=0; ipg<npg; ipg++) {

      Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);
      detJaco = mydetsur(Jaco,normal);
      normal.scale(-1.); // fixme:= This is to compensate a bug in mydetsur
      if (detJaco <= 0.) {
	cout << "bcconv: Jacobian of element " << k << " is negative or null\n"
	     << " Jacobian: " << detJaco << endl ;
	assert(0);
      }

      // WPG is used instead of wpgdet because normal is not a unit vector ; 

      if (comp_mat_res) {

	// state variables and gradient
	p_star = double(tmp1.prod(SHAPE,pcol_star,-1,-1));
	u_star.prod(SHAPE,ucols_star,-1,-1,1);

	// Residual contribution due to integration by parts of pressure gradient
	tmp2.set(normal).scale(p_star);
	tmp3.prod(SHAPE,tmp2,1,2);
	resmom.axpy(tmp3,WPG);

	// Matrix contribution due to integration by parts of pressure gradient
	for (int iloc=1; iloc<=nel; iloc++) {
	  double SHAPE_iloc = SHAPE.get(iloc);
	  for (int jloc=1; jloc<=nel; jloc++) {
	    double SHAPE_jloc = SHAPE.get(jloc);
	    tmp2.set(normal).scale(SHAPE_jloc);
	    tmp2.scale(SHAPE_iloc);
            // matij.set(0.);
	    if (weak_form) {
	      // matij.ir(2,ndof).is(1,1,ndim).axpy(tmp2,SHAPE_iloc);
	      matij.ir(2,ndof).is(1,1,ndim).set(tmp2);
	    } 
	    matij.rs();

	    int il1=(iloc-1)*ndof+1;
	    int il2=il1+ndof-1;
	    int jl1=(jloc-1)*ndof+1;
	    int jl2=jl1+ndof-1;
	    matloc.is(1,il1,il2).is(2,jl1,jl2).axpy(matij,-WPG).rs();
	    
	  }
	}

      } else if (comp_mat) {
	// don't make anything here !!
      } else {
	PetscPrintf(PETSC_COMM_WORLD,
		    "Don't know how to compute jobinfo: %s\n",jobinfo);
	CHKERRQ(ierr);
      }

    }

    // BETO : matloc_prof habria que inicializarlo a cero ???????
    if(comp_mat) {
      matloc_prof.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
    }      

    if (comp_mat_res) {
      if (!weak_form) {
	veccontr.set(0.);
	matloc.set(0.);
      } else {
        veccontr.is(2,1,ndim).set(resmom).rs();
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
