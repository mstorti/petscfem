//__INSERT_LICENSE__
//$Id: bccns_gasflow.cpp,v 1.1.2.1 2005/09/20 00:58:38 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>

#include <src/fastmat2.h>
#include "nsi_tet.h"

extern TextHashTable *GLOBAL_OPTIONS;
#define MAXPROP 100

// ======================================================================================
#define set_state_bcc {						\
	vel_supg.set(vel).rest(v_mesh).rs();			\
        if(axi>0){						\
          vel_supg.setel(0.,axi);				\
        }							\
	u2 = vel_supg.sum_square_all();				\
	velmod = sqrt(u2);					\
	g1=ga-1;						\
	rho_ene = p/g1+0.5*rho*square(velmod);	\
	ene = rho_ene/rho;					\
	entalpy=rho_ene+p;					\
        Cv=Rgas/g1;						\
        int_ene=p/g1/rho;				\
        Cp = ga*Cv;						\
}  

// Advective fluxes							
#define compute_Fa {					\
  flux.set(0.);						\
  flux.ir(1,1).set(vel).scale(rho).rs();		\
  Amom.prod(vel,vel,1,2);				\
  Y.set(0.).add(Amom).scale(rho).axpy(Id_ndim,p);	\
  flux.is(1,vl_indx,vl_indxe).add(Y).rs();		\
  flux.ir(1,ndof).set(vel).scale(entalpy).rs();		\
}

// ======================================================================================
// Adjective Jacobians in primitive basis						      
#define compute_Ajac {								\
  Ajac.set(0.);									\
  Ajac.ir(2,1).ir(3,1).set(vel).rs();						\
  Ajac.ir(2,1).is(3,vl_indx,vl_indxe).set(Id_ndim).scale(rho).rs();		\
  Ajac.ir(2,ndof).ir(3,1).set(vel).scale(0.5*square(velmod)).rs();		\
  Ajac.ir(2,ndof).ir(3,ndof).set(vel).scale(ga/g1).rs();			\
  Ajac.ir(2,ndof).is(3,vl_indx,vl_indxe).prod(vel,vel,1,2).scale(rho).rs();     \
  Ajac.ir(2,ndof).is(3,vl_indx,vl_indxe).axpy(Id_ndim,ga/g1*p+0.5*rho*square(velmod)).rs();   \
  Ajac.is(2,vl_indx,vl_indxe).ir(3,1).prod(vel,vel,1,2).rs();				      \
  Ajac.is(2,vl_indx,vl_indxe).ir(3,ndof).set(Id_ndim).rs();				      \
  Ajac.is(2,vl_indx,vl_indxe).is(3,vl_indx,vl_indxe).prod(vel,Id_ndim,1,2,3).scale(rho).rs(); \
  Y3.prod(vel,Id_ndim,2,1,3).scale(rho).rs();				\
  Ajac.is(2,vl_indx,vl_indxe).is(3,vl_indx,vl_indxe).add(Y3).rs();	\
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "int bcconv_ns_gasflow::ask(char *,int &)"
int bcconv_ns_gasflow::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(comp_mat);
  DONT_SKIP_JOBINFO(comp_res);
  DONT_SKIP_JOBINFO(comp_mat_res);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "bcconv_nsi_tet_asm::assemble"
int bcconv_ns_gasflow::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
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

  int ierr=0, axi;

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

  // Hloc stores the old mesh coordinates
  int nH = nu-ndim;
  FMatrix  Hloc(nel,nH),vloc_mesh(nel,ndim),v_mesh(ndim);

  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  GlobParam *glob_param;
  double *hmin,rec_Dt;
  int ja_hmin;

  // Get arguments from arg_list
  double *locst,*locst2,*retval,*retvalmat;
  if (comp_mat) {
    retvalmat = arg_data_v[0].retval;
  }

#define WAS_SET arg_data_v[ja_hmin].was_set
  if (comp_mat_res) {
    int ja=0;
    locst = arg_data_v[ja++].locst;
    locst2 = arg_data_v[ja++].locst;
    retval = arg_data_v[ja++].retval;
    retvalmat = arg_data_v[ja++].retval;
    hmin = &*(arg_data_v[ja++].vector_assoc)->begin();
    ja_hmin=ja;
    glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
    rec_Dt = 1./glob_param->Dt;
    if (glob_param->steady) rec_Dt=0.;
  }

  // allocate local vecs
  int kdof;
  FMatrix veccontr(nel,ndof),xloc(nel,ndim),locstate(nel,ndof),
    locstate2(nel,ndof),xpg;

  nen = nel*ndof;

  FMatrix matloc(nen,nen),matloc_prof(nen,nen);

  // Physical properties
  int iprop=0, elprpsindx[MAXPROP]; double propel[MAXPROP];

  //o Add axisymmetric version for this particular elemset.
  TGETOPTDEF_S(thash,string,axisymmetric,none);
  assert(axisymmetric.length()>0);
  if (axisymmetric=="none") axi=0;
  else if (axisymmetric=="x") axi=1;
  else if (axisymmetric=="y") axi=2;
  else if (axisymmetric=="z") axi=3;
  else {
    PetscPrintf(PETSC_COMM_WORLD,
		"Invalid value for \"axisymmetric\" option\n"
		"axisymmetric=\"%s\"\n",axisymmetric.c_str());
    PetscFinalize();
    exit(0);
  }

  // Trapezoidal method parameter.
  TGETOPTDEF(GLOBAL_OPTIONS,double,alpha,1.);

  //o ALE_flag : flag to ON ALE computation
  SGETOPTDEF(int,ALE_flag,0);
  //o indx_ALE_xold : pointer to old coordinates in
  //  NODEDATA array excluding the first "ndim" values
  SGETOPTDEF(int,indx_ALE_xold,1);

  // gamma coeficient
  SGETOPTDEF(double,ga,1.4);
  // constant of a particular gas for ideal gas law (state equation for the gas)
  SGETOPTDEF(double,Rgas,287.);

  double pi = 4*atan(1.0);

  int vl_indx=2, vl_indxe = vl_indx+ndim-1;

  int nprops=iprop;

  FastMat2 rhocol,rhocol_new,rhocol_star;
  FastMat2 pcol,pcol_new,pcol_star;
  FastMat2 vel,ucols,ucols_new,ucols_star,u_star,vel_supg;
  FastMat2 tmp1,tmp2,tmp3,tmp4,tmp5;

  FastMat2 Amom,Y(2,ndim,ndim),Y3,Ajac(3,ndim,ndof,ndof),Ajac_n(2,ndof,ndof),fluxn(1,ndof),
    flux(2,ndof,ndim),Id_ndim(2,ndim,ndim);

  double velmod,u2;
  double g1,rho_ene,ene,entalpy,Cv,int_ene,Cp;
  double rho,p; 

  Id_ndim.eye();

  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  GPdata gp_data(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);

  // Definiciones para descargar el lazo interno
  double detJaco,rho_star,p_star,wpgdet;

  int elem, ipg,node, jdim, kloc,lloc,ldof;

  FMatrix Jaco(ndimel,ndim),normal(ndim);

  matloc_prof.set(1.);
   
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
      if(nH>0) Hloc.ir(1,kloc+1).set(&NODEDATA(node-1,0)+ndim);
    }
    xloc.rs();
    Hloc.rs();

    // tenemos el estado locstate2 <- u^n
    //                   locstate  <- u^*
    if (comp_mat_res) {
      locstate.set(&(LOCST(ielh,0,0)));
      locstate2.set(&(LOCST2(ielh,0,0)));
    }

    matloc.set(0.);
    veccontr.set(0.);

    ucols.set(locstate2.is(2,vl_indx,vl_indxe));
    pcol.set(locstate2.rs().ir(2,ndof));
    rhocol.set(locstate2.rs().ir(2,1));
    locstate2.rs();

    ucols_new.set(locstate.is(2,vl_indx,vl_indxe));
    pcol_new.set(locstate.rs().ir(2,ndof));
    rhocol_new.set(locstate.rs().ir(2,1));
    locstate.rs();
      
    ucols_star.set(ucols_new).scale(alpha).axpy(ucols,1-alpha);
    pcol_star.set(pcol_new).scale(alpha).axpy(pcol,1-alpha);
    rhocol_star.set(rhocol_new).scale(alpha).axpy(rhocol,1-alpha);

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    // nodal computation of mesh velocity
    if (ALE_flag) {
      assert(nH >= ndim);
      assert(indx_ALE_xold >= nH+1-ndim);
      Hloc.is(2,indx_ALE_xold,indx_ALE_xold+ndim-1);
      vloc_mesh.set(xloc).rest(Hloc).scale(rec_Dt).rs();
      Hloc.rs();
    }

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
      
      if (comp_mat_res) {
	
	// current state
	rho_star = double(tmp4.prod(SHAPE,rhocol_star,-1,-1));
	p_star = double(tmp4.prod(SHAPE,pcol_star,-1,-1));
	u_star.prod(SHAPE,ucols_star,-1,-1,1);

	v_mesh.set(0.0);
	if (ALE_flag) {
	  v_mesh.prod(SHAPE,vloc_mesh,-1,-1,1);
	}
	
	// set state to current state
	rho=rho_star;p=p_star;vel.set(u_star);
	set_state_bcc;

	compute_Fa;
	
	fluxn.prod(flux,normal,1,-1,-1);
	tmp1.prod(SHAPE,fluxn,1,2);
	veccontr.axpy(tmp1,WPG);
	
	compute_Ajac;
	Ajac_n.prod(Ajac,normal,-1,1,2,-1);
	
	//	tmp2.prod(SHAPE,Ajac_n,1,2,3); // (nel,ndof,ndof)
	//	tmp3.prod(tmp2,SHAPE,1,2,4,3); // (nel,ndof,nel,ndof)

	//        tmp5.set(tmp3.resize(2,nen,nen));
	//	matloc.axpy(tmp5,-WPG);	       // (nen,nen)	

	for (int iloc=1; iloc<=nel; iloc++) {
	  double SHAPE_iloc = SHAPE.get(iloc);
	  for (int jloc=1; jloc<=nel; jloc++) {
	    double SHAPE_jloc = SHAPE.get(jloc);
	    int il1=(iloc-1)*ndof+1;
	    int il2=iloc*ndof;
	    int jl1=(jloc-1)*ndof+1;
	    int jl2=jloc*ndof;
	    double vaux = SHAPE_iloc*SHAPE_jloc;
	    matloc.is(1,il1,il2).is(2,jl1,jl2).axpy(Ajac_n,-vaux*WPG).rs();	    
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

    if(comp_mat) {
      matloc_prof.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
    }

    if (comp_mat_res) {
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
