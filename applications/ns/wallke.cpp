//__INSERT_LICENSE__
//$Id: wallke.cpp,v 1.11 2001/06/29 16:34:24 mstorti Exp $
#include "../../src/fem.h"
#include "../../src/utils.h"
#include "../../src/readmesh.h"
#include "../../src/getprop.h"

#include "../../src/fastmat2.h"
#include "nsi_tet.h"

extern TextHashTable *GLOBAL_OPTIONS;
extern int MY_RANK,SIZE;
extern int TSTEP; //debug:=
#define MAXPROP 100

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retvalmat,iele,j,nel,k,ndof,p,nel,q,ndof)
#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define IDENT(j,k) (ident[ndof*(j)+(k)]) 
#define JDOFLOC(j,k) VEC2(jdofloc,j,k,ndof)

// This is required!! See `Thinking in C++, 2nd ed. Volume 1', 
// (http://www.MindView.net)
// Chapter 15: `Polymorphism &  Virtual Functions', 
// paragraph `Pure virtual destructors' 
WallFun::~WallFun() {}; 

void WallFun1::w(double yp,double &f,double &fprime) {
  if (yp<0) {
    PETSCFEM_ERROR0("y+<0. Invalid value for y+\n");
  } if (yp<5) {
    f = yp;
    fprime = 1.;
  } else if (yp<30) {
    f = 5.0*log(yp)-3.05;
    fprime = 5.0/yp;
  } else {
    f = 2.5*log(yp)+5.5;
    fprime = 2.5/yp;
//      f = 1./Chi * log(E_star*yp);
//      fprime = 1./Chi/yp;
  }
}

#undef __FUNC__
#define __FUNC__ "int NonLinearRes::ask(const char *jobinfo,int &skip_elemset)"
int wallke::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(comp_mat);
  DONT_SKIP_JOBINFO(comp_mat_res);
  return 0;
}
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "wall::assemble"
int wallke::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
			  Dofmap *dofmap,const char *jobinfo,int myrank,
			  int el_start,int el_last,int iter_mode,
			  const TimeData *time_data) {

  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_mat_res);

  int ierr=0;

  int npg,ndim;
  ierr = get_int(thash,"npg",&npg); CHKERRA(ierr);
  ierr = get_int(thash,"ndim",&ndim); CHKERRA(ierr);
  double i_nel = 1./double(nel);

  int ndimel = ndim-1;
  int nen = nel*ndof;

  int nu=nodedata->nu;

  // Physical properties
  int iprop=0, elprpsindx[MAXPROP]; double propel[MAXPROP];
  DEFPROPN(u_wall,ndim);
  int nprops=iprop;

  // Get arguments from arg_list
  double *locst,*locst2,*retval,*retvalmat;
  Elemset *elemset;

  //o The $y^+$ coordinate of the computational boundary
  SGETOPTDEF(double,y_wall_plus,30.);
  //o The $y$ (normal) coordinate of the computational boundary. 
  // Only for laminar computations.
  SGETOPTDEF(double,y_wall,1e-3);
  //o Mask for using laminar relation (\verb+turbulence_coef=0+). 
  SGETOPTDEF(double,turbulence_coef,1.);
  //o Use lumped mass matric for the wall element contribution. Avoids
  // oscillations due to ``reactive type'' wall contributions. 
  SGETOPTDEF(int,lumped_wallke,1);

  SGETOPTDEF(double,viscosity,0.); //o
  PETSCFEM_ASSERT0(viscosity>0.,
		   "Not entered `viscosity' property or "
		   "null value  entered.\n");
  WallFun1 wf1(this);
  WallFun *wf = &wf1;
  wf->init();

  double fwall,fprime;
  if (comp_mat_res) {
    wf->w(y_wall_plus,fwall,fprime);

    locst = arg_data_v[0].locst;
    locst2 = arg_data_v[1].locst;
    if (comp_mat_res) {
      retval = arg_data_v[2].retval;
      retvalmat = arg_data_v[3].retval;
    }
  }

  if (comp_mat) {
    retvalmat = arg_data_v[0].retval;
  }

  // allocate local vecs
  int kdof;
  FMatrix veccontr(nel,ndof),xloc(nel,ndim),xc(ndim),locstate(nel,ndof),
    locstate2(nel,ndof),xpg; 

  nen = nel*ndof;
  FastMat2 matloc(4,nel,ndof,nel,ndof),lmass(1,nel ),u_wall(1,ndim);
  FMatrix matlocmom(nel,nel);

  double rho=1.;
  // Trapezoidal method parameter. 
  TGETOPTDEF(GLOBAL_OPTIONS,double,alpha,1.);    //nd

  // Gauss Point data
  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  
  GPdata gp_data(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);

  // Definiciones para descargar el lazo interno
  double detJaco,p_star,wpgdet;

  int elem, ipg,node, jdim, kloc,lloc,ldof;
    
  FMatrix Jaco(ndimel,ndim),resmom(nel,ndim),normal(ndim);
          
  FMatrix grad_p_star(ndim),u,u_star,ucols,ucols_new,ucols_star,
    pcol_star,pcol_new,pcol;

  FMatrix matloc_prof(nen,nen),uc(ndim),tmp1,tmp2,tmp3,tmp4,
    tmp5,seed(ndof,ndof);

//    if (comp_mat_res) {
//      seed.eye().setel(0.,ndof,ndof);
//    }
  seed.set(0.);
  for (int j=1; j<=ndim; j++) 
    seed.setel(1.,j,j);

  if (comp_mat) {
    matloc_prof.set(1.);
  }

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  int ielh=-1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    FastMat2::reset_cache();

    ielh++;
    load_props(propel,elprpsindx,nprops,&(ELEMPROPS(k,0)));
    u_wall.set(propel+u_wall_indx);

    // Load local node coordinates in local vector
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(k,kloc);
      xloc.ir(1,kloc+1).set(&NODEDATA(node-1,0));
    }
    xloc.rs();

    // locstate2 <- u^n
    // locstate  <- u^*
    if (comp_mat_res ) {
      locstate.set(&(LOCST(ielh,0,0)));
      locstate2.set(&(LOCST2(ielh,0,0)));
    }

    matlocmom.set(0.);
    matloc.set(0.);
    veccontr.set(0.);

    ucols.set(locstate2.is(2,1,ndim));
    locstate2.rs();

    ucols_new.set(locstate.is(2,1,ndim));
    locstate.rs();

    ucols_star.set(ucols_new).scale(alpha).axpy(ucols,1-alpha);

    if(comp_mat) {
      matloc_prof.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
      continue;
    }      

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])
    
    if (comp_mat_res) {
      veccontr.is(2,1,ndim);
      matloc.is(2,1,ndim).is(4,1,ndim);
      lmass.set(0.)       ;

      // loop over Gauss points
      // Guarda que hay que invertir la direccion de la traccion!!!!
      for (ipg=0; ipg<npg; ipg++) {

	Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);
	detJaco = mydetsur(Jaco,normal);
	wpgdet = detJaco * WPG;
	normal.scale(-1.); // fixme:= This is to compensate a bug in mydetsur

	u_star.prod(SHAPE,ucols_star,-1,-1,1).rest(u_wall);
	double Ustar = sqrt(u_star.sum_square_all());
	double gfun,gprime;
	gprime = rho / (fwall*fwall);
	gfun = gprime * Ustar;
	if (turbulence_coef == 0.) {
	  gprime = 0.;
	  gfun = rho * viscosity / y_wall;
	}

	tmp1.prod(SHAPE,SHAPE,1,2);
	matlocmom.axpy(tmp1,gfun*wpgdet);
	if (lumped_wallke) 
	  lmass.axpy(SHAPE,wpgdet);

	tmp2.prod(SHAPE,u_star,1,2);
	tmp3.set(tmp2).scale(wpgdet*gprime/Ustar);
	tmp4.prod(tmp2,tmp3,1,2,3,4);
	matloc.add(tmp4);

	veccontr.axpy(tmp2,-gfun*wpgdet);
      }

      if (lumped_wallke) {
	matloc.reshape(2,nel*ndof,nel*ndof);
	matloc.reshape(4,nel,ndof,nel,ndof);

	matlocmom.set(0.);
	matloc.rs().set(0.).is(2,1,ndim).is(4,1,ndim);
	veccontr.rs().set(0.).is(2,1,ndim);

	for (int j=1; j<=nel; j++) {
	  ucols_star.ir(1,j);
	  double Ustar = sqrt(ucols_star.sum_square_all());
	  double gfun,gprime, Omega_j;
	  Omega_j = lmass.get(j);
	  gprime = rho / (fwall*fwall);
	  gfun = gprime * Ustar;
	  if (turbulence_coef == 0.) {
	    gprime = 0.;
	    gfun = rho * viscosity / y_wall;
	  }

	  matlocmom.setel(Omega_j*gfun,j,j);

	  matloc.ir(1,j).ir(3,j)
	    .prod(u_star,u_star,1,2)
	    .scale(wpgdet*gprime/Ustar);

	  veccontr.ir(1,j).set(ucols_star).scale(-gfun*Omega_j);
	}
      }

      matlocmom.rs();
      tmp5.prod(matlocmom,seed,1,3,2,4);

      veccontr.rs();
      matloc.rs().add(tmp5);

      matloc.reshape(2,nel*ndof,nel*ndof);
      matloc.reshape(4,nel,ndof,nel,ndof);

      veccontr.export_vals(&(RETVAL(ielh,0,0)));
      matloc.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
    }

  }

  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  return 0;
}
#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
