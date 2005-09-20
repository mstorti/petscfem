//__INSERT_LICENSE__
//$Id: bccnsitetasm.cpp,v 1.1.2.1 2005/09/20 00:58:38 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>

#include <src/fastmat2.h>
#include "nsi_tet.h"
#include "nsifunaux.h"

extern TextHashTable *GLOBAL_OPTIONS;
#define MAXPROP 100


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "int bcconv_nsi_tet_asm::ask(char *,int &)"
int bcconv_nsi_tet_asm::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(comp_mat);
  DONT_SKIP_JOBINFO(comp_res);
  DONT_SKIP_JOBINFO(comp_mat_res);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "bcconv_nsi_tet_asm::assemble"
int bcconv_nsi_tet_asm::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
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
  int npg,ndim,nphases;
  ierr = get_int(thash,"npg",&npg); CHKERRA(ierr);
  ierr = get_int(thash,"ndim",&ndim); CHKERRA(ierr);
  ierr = get_int(thash,"nphases",&nphases); CHKERRA(ierr);

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

  //o Solve coupled
  SGETOPTDEF(int,coupled,0);

  // allocate local vecs
  int kdof;
  FMatrix veccontr(nel,ndof),xloc(nel,ndim),locstate(nel,ndof),
    locstate2(nel,ndof),xpg;

  if (ndof != ndim+1+nphases) {
    PetscPrintf(PETSC_COMM_WORLD,"ndof != ndim+1+nphases\n"); CHKERRA(1);
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
  double rho_m,arho_l,arho_g,alpha_l,alpha_g,vslip_m;

  FastMat2 arho_g_vp,v_g_vp,vslip_user_vp,vslip_vp,vslip_m_vp,rho_g_vp;
  FastMat2 vfcols,vfcols_new,vfcols_star,Id_vp(2,nphases,nphases);
  FastMat2 vf,vf_star,alpha_g_vp;

  FastMat2 d_bubble_vp,tmp2_cont, tmp3_cont,tmp_g,tmp2_g,tmp3_g;
  double vslip, rho_g;

  // Gauss Point data
  //    char *geom;
  //    thash->get_entry("geometry",geom);

  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  GPdata gp_data(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);

  // Definiciones para descargar el lazo interno
  double detJaco,p_star,wpgdet;

  int elem, ipg,node, jdim, kloc,lloc,ldof;

  FMatrix Jaco(ndimel,ndim),resmom(nel,ndim),normal(ndim),matij(ndim+1,ndim+1);
  FMatrix rescont(nel),res_alpha_g(nel,nphases);

  FMatrix grad_p_star(ndim),u,u_star,ucols,ucols_new,ucols_star,
    pcol_star,pcol_new,pcol,tmp1,tmp2,tmp3,tmp4;

  FastMat2 bcconv_factor(1,ndof);

  TGETOPTDEF(GLOBAL_OPTIONS,double,rho_l,0.);
  assert(rho_l>0);

  TGETOPTDEF(GLOBAL_OPTIONS,int,use_modified_slag_vslip,0);

  //  phases density vector
  rho_g_vp.resize(1,nphases);
  // rho_g_vp.set(rho_g);
  ierr = get_double(GLOBAL_OPTIONS,"rho_phases",rho_g_vp.storage_begin(),1,nphases);

  //  Bubble diameter vector
  d_bubble_vp.resize(1,nphases);
  d_bubble_vp.set(0.0);
  ierr = get_double(GLOBAL_OPTIONS,"d_bubble_phases",
		    d_bubble_vp.storage_begin(),1,nphases);

  //  slip velocity vector
  vslip_user_vp.resize(1,nphases);
  vslip_user_vp.set(0.0);
  ierr = get_double(GLOBAL_OPTIONS,"vslip_user_phases",
		    vslip_user_vp.storage_begin(),1,nphases);

  //o Direction of gravity
  TGETOPTDEF(thash,int,g_dir,ndim);
  
  bcconv_factor.set(1.);

  //  ierr = get_double(GLOBAL_OPTIONS,"bcconv_factor",
  //	    bcconv_factor.storage_begin(),1,ndof);  

  ierr = get_double(thash,"bcconv_factor",
		    bcconv_factor.storage_begin(),1,ndof);  

  v_g_vp.resize(2,ndim,nphases);
  alpha_g_vp.resize(1,nphases);
  vslip_vp.resize(1,nphases);

  FMatrix matloc_prof(nen,nen),seed,one_nel;

  //  if (comp_mat) matloc_prof.set(1.);

  // ====================================

    if (coupled) {
      matloc_prof.set(1.);
    } else {      
      
      seed.resize(2,ndof,ndof);
      seed.set(0.);
      one_nel.resize(2,nel,nel);
      one_nel.set(0.);
      matloc_prof.resize(2,nen,nen);
      
#ifndef ADD_GRAD_DIV_U_TERM
      for (int jj=1; jj<=ndim; jj++) {
	seed.setel(1.,jj,jj);
	seed.setel(1.,jj,ndim+1);
	seed.setel(1.,ndim+1,jj);
      }
      seed.setel(1.,ndim+1,ndim+1);
#else
      seed.is(1,1,ndim+1).is(2,1,ndim+1).set(1.0).rs();
#endif

      // disperse phases
      for (int j=ndim+2; j<=ndim+1+nphases; j++) {
	seed.setel(1.,j,j);
      }
      one_nel.set(1.);
      matloc_prof.kron(one_nel,seed);
    }

    // ====================================


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
    if (comp_mat_res) {
      locstate.set(&(LOCST(ielh,0,0)));
      locstate2.set(&(LOCST2(ielh,0,0)));
    }

    matlocmom.set(0.);
    matloc.set(0.);
    veccontr.set(0.);
    resmom.set(0.);
    rescont.set(0.);
    res_alpha_g.set(0.);

    ucols.set(locstate2.is(2,1,ndim));
    pcol.set(locstate2.rs().ir(2,ndim+1));
    vfcols.set(locstate2.rs().is(2,ndim+2,ndim+1+nphases));
    locstate2.rs();

    ucols_new.set(locstate.is(2,1,ndim));
    pcol_new.set(locstate.rs().ir(2,ndim+1));
    vfcols_new.set(locstate.rs().is(2,ndim+2,ndim+1+nphases));
    locstate.rs();

    ucols_star.set(ucols_new).scale(alpha).axpy(ucols,1-alpha);
    pcol_star.set(pcol_new).scale(alpha).axpy(pcol,1-alpha);
    vfcols_star.set(vfcols_new).scale(alpha).axpy(vfcols,1-alpha);

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    matij.set(0.);

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

	// state variables and gradient
	p_star = double(tmp1.prod(SHAPE,pcol_star,-1,-1));
	u_star.prod(SHAPE,ucols_star,-1,-1,1);
	vf_star.prod(SHAPE,vfcols_star,-1,-1,1);

	alpha_g_vp.set(vf_star);

	rho_m = compute_rho_m(rho_g_vp,arho_g_vp,alpha_g_vp,alpha_l,rho_l,nphases);

	compute_vel_g(u_star,vslip_vp,vslip_user_vp,vslip_m_vp,v_g_vp,rho_g_vp,
		      rho_m,g_dir,d_bubble_vp,nphases,use_modified_slag_vslip,
		      vf_star);
	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

	// Residual contribution due to integration by parts of div(rho_m*Vm)
	tmp2_cont.prod(normal,u_star,-1,-1).scale(rho_m);
	tmp3_cont.set(SHAPE).scale(double(tmp2_cont));
	rescont.axpy(tmp3_cont,-WPG);

	// Residual contribution due to integration by parts of pressure gradient
	tmp2.set(normal).scale(p_star);
	tmp3.prod(SHAPE,tmp2,1,2);
	resmom.axpy(tmp3,WPG);

	// Residual cont due to integration by parts of div(alpha_k*V_k)

	for (int j=1; j<=nphases; j++) {
	  v_g_vp.ir(2,j);
	  res_alpha_g.ir(2,j);
	  tmp2_g.prod(normal,v_g_vp,-1,-1).scale(double(vf_star.get(j)));
	  tmp3_g.set(SHAPE).scale(double(tmp2_g));
	  res_alpha_g.axpy(tmp3_g,WPG);
	}
	res_alpha_g.rs();
	v_g_vp.rs();

	resmom.ir(2,abs(g_dir));
	for (int j=1; j<=nphases; j++) {
	  rho_g = rho_g_vp.get(j);
	  alpha_g = alpha_g_vp.get(j);
	  vslip_m = vslip_m_vp.get(j);
	  double vaux = rho_g*vslip_m*vslip_m;
	  double no_dot_g = normal.get(abs(g_dir));
	  resmom.axpy(SHAPE,WPG*alpha_g*vaux*no_dot_g);
	}
	resmom.rs();

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

	    // contribucion a la matriz de la debilitacion de la ec.
	    // de continuidad Ni*div(rho_m*V_m)
	    matloc.ir(1,il2).is(2,jl1,jl2-1).axpy(tmp2,rho_m*WPG).rs();

	    // contribucion a las ecuaciones de la fase dispersa
	    // por la debilitacion del termino convectivo
	    for (int j=1; j<=nphases; j++) {
	      v_g_vp.ir(2,j);
	      double NiNjvn = double(tmp_g.prod(v_g_vp,tmp2,-1,-1));
	      il1 = (iloc-1)*ndof+ndim+1+j;
	      jl1 = (jloc-1)*ndof+ndim+1+j;
	      matloc.addel(-WPG*NiNjvn,il1,jl1);	      
	      if (coupled) {
		// contribucion del termino que surje al debilitar la ecuacion 
		// de continuidad porque rho_m depende de alpha_g
		rho_g = rho_g_vp.get(j);
		il1 = (iloc-1)*ndof+ndim+1;
		jl1 = (jloc-1)*ndof+ndim+1+j;
		matloc.addel(WPG*(rho_g-rho_l)*NiNjvn,il1,jl1);
	      }
	    }

	    if (coupled) {	      
	      // contribucion a las ecuaciones de momento de la mezcla
	      // por la debilitacion del termino extra que aparece al
	      // expresar la ecuacion de momento de la mezcla en funcion
	      // de la velocidad de la mezcla 
	      // div(sum(rho_k*alpha_k*dyad(vslip_k,vslip_k)))
	      for (int j=1; j<=nphases; j++) {
		rho_g = rho_g_vp.get(j);
		vslip_m = vslip_m_vp.get(j);
		double vaux = rho_g*vslip_m*vslip_m;
		double NiNjno_dot_g = tmp2.get(abs(g_dir));
		il1 = (iloc-1)*ndof+abs(g_dir);
		jl1 = (jloc-1)*ndof+ndim+1+j;
		matloc.addel(-WPG*vaux*NiNjno_dot_g,il1,jl1);
	      }
	    }
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
        veccontr.ir(2,1).set(rescont).rs();
        veccontr.is(2,ndim+2,ndim+1+nphases).set(res_alpha_g).rs();

	/*
	printf(" DEBUG --- > bcconv \n");
	xloc.print("xloc :");
	res_alpha_g.print("res_alpha_g :");
	*/

      }

    // apply bcconv_factor to mask residual and matrix contributions
    for (int j=1; j<=ndof; j++) {      
      veccontr.ir(2,j).scale(bcconv_factor.get(j)).rs();
      matloc.is(1,j,nen,ndof).scale(bcconv_factor.get(j)).rs();
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
