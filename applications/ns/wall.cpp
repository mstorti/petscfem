//__INSERT_LICENSE__
//$Id: wall.cpp,v 1.20 2004/10/01 00:55:22 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>

#include <src/fastmat2.h>
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
#define ELEMPROPS_ADD(j,k) VEC2(elemprops_add,j,k,nelprops_add)
#define SHEAR_VEL(j) ELEMPROPS_ADD(j,0)
#define IDENT(j,k) (ident[ndof*(j)+(k)])
#define JDOFLOC(j,k) VEC2(jdofloc,j,k,ndof)


extern vector<double> data_pts;
extern vector<ElemToPtr> elemset_pointer;


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "int wall::ask(char *,int &)"
int wall::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  int ierr;

  // Use LES or not.
  TGETOPTDEF(GLOBAL_OPTIONS,int,LES,0); //nd
  if (!LES) return 0;

  DONT_SKIP_JOBINFO(build_nneighbor_tree);
  DONT_SKIP_JOBINFO(comp_shear_vel);
  DONT_SKIP_JOBINFO(comp_mat_res);
  DONT_SKIP_JOBINFO(comp_mat);

  // Initialize to 0. all shear_velocities
  if (!strcmp(jobinfo,"comp_shear_vel")) {
    for (int k=0; k<nelem; k++) SHEAR_VEL(k)=0.;
  }

  if (!strcmp(jobinfo,"communicate_shear_vel") && SIZE>1 ) {
    // I dont't know if I can do the reduce `in place',
    // meanwhile, I use a temporary storage
    double *recv_buff = new double[nelem];
    ierr = MPI_Allreduce((void *)elemprops_add,
			 (void *)recv_buff,nelem,MPI_DOUBLE,
			 MPI_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);
    for (int j=0; j<nelem; j++) elemprops_add[j] = recv_buff[j];
    delete[] recv_buff;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void wall::initialize() {
  int ierr;
  ierr = get_int(thash,"ndim",&ndim);
  assert(!ierr);
  // convert pointers
  data_pts.resize(ndim*nelem); // Not implemented yet
				// Number of wall elemsets >1
  elemset_pointer.push_back(ElemToPtr(nelem,this));
  assert(elemset_pointer.size()==1); // Not implemented yet
				// Number of wall elemsets >1
				// see above!!
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ ""
void wall_fun(double yp,double &f,double &fprime) {

  static double a1=5.0, b1=-3.05, a2=2.5, b2=5.5,
    yp1=5, yp2=30;

  if (yp < yp1) {
    f = yp;
    fprime = 1.;
  } else if (yp < yp2) {
    f = (a1*log(yp)+b1)/yp;
    fprime = (a1 - b1 - a1*log(yp))/(yp*yp);
  } else {
    f = (a2*log(yp)+b2)/yp;
    fprime = (a2 - b2 - a2*log(yp))/(yp*yp);
  }
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "wall::assemble"
int wall::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
			  Dofmap *dofmap,const char *jobinfo,int myrank,
			  int el_start,int el_last,int iter_mode,
			  const TimeData *time_data) {

  // wait_from_console("entra a wall::assemble");
  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_mat_res);
  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(build_nneighbor_tree);
  GET_JOBINFO_FLAG(comp_shear_vel);

  int ierr=0;

  int npg;
  ierr = get_int(thash,"npg",&npg); CHKERRA(ierr);
  double i_nel = 1./double(nel);

  int ndimel = ndim-1;
  int nen = nel*ndof;

  int nu=nodedata->nu;

  // Get arguments from arg_list
  double *locst,*locst2,*retval,*retvalmat;

  //o The $y^+$ coordinate of the computational boundary
  TGETOPTDEF(thash,double,y_wall_plus,25.);
  double fwall,fprime;
  if (comp_shear_vel || comp_mat_res) {
    wall_fun(y_wall_plus,fwall,fprime);

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
  FastMat2 matloc(4,nel,ndof,nel,ndof);
  FMatrix matlocmom(nel,nel);

  //  double rho=1.;
  //o Density
  TGETOPTDEF(thash,double,rho,1.);

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

  FMatrix matloc_prof(nen,nen),uc(ndim),tmp1,tmp2,tmp3,tmp4,tmp5,seed(ndof,ndof);

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
    // Load local node coordinates in local vector
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(k,kloc);
      xloc.ir(1,kloc+1).set(&NODEDATA(node-1,0));
    }
    xloc.rs();

    // tenemos el estado locstate2 <- u^n
    //                   locstate  <- u^*
    if (comp_mat_res || comp_shear_vel) {
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

    if (comp_shear_vel) {
      double vel = sqrt(uc.sum(ucols_star,-1,1).sum_square_all())*i_nel;
      SHEAR_VEL(k) = vel/fwall;
      continue;
    }

    if(comp_mat) {
      matloc_prof.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
      continue;
    }

    if (build_nneighbor_tree) {
      // This doesn't compile under RH 5.2
      xc.sum(xloc,-1,1).scale(i_nel);
      double *xc_ = xc.storage_begin();
      for (int j=0; j<ndim; j++) data_pts[k*ndim+j] = xc_[j];
      continue;
    }

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    if (comp_mat_res) {
      // loop over Gauss points
      veccontr.is(2,1,ndim);
      matloc.is(2,1,ndim).is(4,1,ndim);

      // Guarda que hay que invertir la direccion de la traccion!!!!
      for (ipg=0; ipg<npg; ipg++) {

	Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);
	detJaco = Jaco.detsur(&normal);
	normal.scale(-1.); // fixme:= This is to compensate a bug in mydetsur

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
