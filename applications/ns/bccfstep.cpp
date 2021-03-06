//__INSERT_LICENSE__
//$Id merge-with-petsc-233-55-g52bd457 Fri Oct 26 13:57:07 2007 -0300$
  
#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>

#include <src/fastmat2.h>
#include "nsi_tet.h"
#include "fracstep.h"

extern TextHashTable *GLOBAL_OPTIONS;
#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int bcconv_fstep_fm2::ask(char *,int &)"
int bcconv_fstep_fm2::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(comp_res_poi);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "bcconv_fstep_fm2::assemble"
int bcconv_fstep_fm2::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
			       Dofmap *dofmap,const char *jobinfo,int myrank,
			       int el_start,int el_last,int iter_mode,
			       const TimeData *time_data) {

  GET_JOBINFO_FLAG(comp_res_poi);

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

  //o Reverse sense of normal. Useful if, for instance,
  //  you generated, by mistake, the wrong sense of
  //  numeration of nodes in the connectivities. 
  TGETOPTDEF(thash,int,reverse_normal,0);

  int ndimel = ndim-1;
  int nen = nel*ndof;

  // Unpack Dofmap
  int *ident,nnod;
  nnod = dofmap->nnod;

  // Unpack nodedata
  int nu=nodedata->nu;
  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  // Get arguments from arg_list
  double *locst=NULL,*locst2=NULL,*retval=NULL,*retvalmat=NULL;

  GlobParam *glob_param=NULL;
  double Dt=NAN;

  if (comp_res_poi) {
    int ja=0;
    locst = arg_data_v[ja++].locst;
    locst2 = arg_data_v[ja++].locst;
    retval = arg_data_v[ja++].retval;
    glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
    Dt = glob_param->Dt;
  }

  //o Use a weak form for the gradient of pressure term. 
  SGETOPTDEF(int,weak_form,1);

  // allocate local vecs
  int kdof;
  FMatrix veccontr(nel,ndof),xloc(nel,ndim),locstate(nel,ndof),
    locstate2(nel,ndof),xpg; 

#if 0				// allow nsikeps elemset
  if (ndof != ndim+1) {
    PetscPrintf(PETSCFEM_COMM_WORLD,"ndof != ndim+1\n"); CHKERRA(1);
  }
#endif

  nen = nel*ndof;
  //  FMatrix matloc(nen,nen), matlocmom(nel,nel);

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
  GPdata gp_data(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);

  // Definiciones para descargar el lazo interno
  double detJaco,p_star,wpgdet;

  int elem=0, ipg,node, jdim, kloc,lloc,ldof;
    
  FMatrix Jaco(ndimel,ndim),resmom(nel,ndim),normal(ndim),matij(ndim+1,ndim+1);
  FMatrix rescont(nel);
          
  FMatrix grad_p_star(ndim),u,u_star,ucols,ucols_new,ucols_star,
    pcol_star,pcol_new,pcol,tmp1,tmp2,tmp3,tmp4;

  FMatrix matloc_prof(nen,nen);

  //  if (comp_prof) matloc_prof.set(1.);
    
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
    if (comp_res_poi) {
      locstate.set(&(LOCST(ielh,0,0)));
      locstate2.set(&(LOCST2(ielh,0,0)));
    }

    //    matlocmom.set(0.);
    //    matloc.set(0.);
    veccontr.set(0.);
    rescont.set(0.);

    ucols.set(locstate2.is(2,1,ndim));
    //    pcol.set(locstate2.rs().ir(2,ndim+1));
    locstate2.rs();

    ucols_new.set(locstate.is(2,1,ndim));
    //    pcol_new.set(locstate.rs().ir(2,ndim+1));
    locstate.rs();

    ucols_star.set(ucols_new).scale(alpha).axpy(ucols,1-alpha);
    //    pcol_star.set(pcol_new).scale(alpha).axpy(pcol,1-alpha);

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    //    matij.set(0.);

    // loop over Gauss points
    for (ipg=0; ipg<npg; ipg++) {

      Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);
      detJaco = Jaco.detsur(&normal);
      normal.scale(-1.); // fixme:= This is to compensate a bug in mydetsur
      if (detJaco<=0.) {
	detj_error(detJaco,elem);
	set_error(1);
      }

      // WPG is used instead of wpgdet because normal is not a unit vector ; 

      if (comp_res_poi) {

	// state variables and gradient
	//	p_star = double(tmp1.prod(SHAPE,pcol_star,-1,-1));
	u_star.prod(SHAPE,ucols_star,-1,-1,1);

	if (reverse_normal) normal.scale(-1.);

	// Residual contribution due to integration by parts of pressure gradient
	//	tmp2.set(normal).scale(p_star);
        tmp2.prod(u_star,normal,-1,-1);
	//	tmp3.prod(SHAPE,tmp2,1,2);
	tmp3.set(SHAPE).scale(-tmp2*rho/Dt);
	rescont.axpy(tmp3,WPG);

	/*
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
	      matij.ir(2,ndim+1).is(1,1,ndim).set(tmp2);
	    } 
	    matij.rs();

	    int il1=(iloc-1)*ndof+1;
	    int il2=il1+ndim;
	    int jl1=(jloc-1)*ndof+1;
	    int jl2=jl1+ndim;
	    matloc.is(1,il1,il2).is(2,jl1,jl2).axpy(matij,-WPG).rs();
	    
	  }
	}
	*/
      } else {
	PetscPrintf(PETSCFEM_COMM_WORLD,
		    "Don't know how to compute jobinfo: %s\n",jobinfo);
	CHKERRQ(ierr);
      }

    }

    /*
    if(comp_prof) {
      matloc_prof.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
    }      
    */

    if (comp_res_poi) {
      if (!weak_form) {
	veccontr.set(0.);
	//	matloc.set(0.);
      } else {
        veccontr.ir(2,ndof).set(rescont).rs();
      }
	
      veccontr.export_vals(&(RETVAL(ielh,0,0)));
      //      matloc.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
    }
  }
      
  FastMat2::void_cache();
  FastMat2::deactivate_cache();

  return 0;
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
