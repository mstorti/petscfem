//__INSERT_LICENSE__
//$Id: genload.cpp,v 1.1 2001/05/21 17:44:18 mstorti Exp $
 
#include "../../src/fem.h"
#include "../../src/utils.h"
#include "../../src/readmesh.h"
//#include "fastmat.h"
#include "newmat.h"
#include "../../src/getprop.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int advective::ask(char *,int &)"
int GenLoad::ask(const char *jobinfo,int &skip_elemset) {
   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_prof);
   return 0;
}

NSGETOPTDEF_ND(hfilm,);


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "GenLoad::new_assemble"
void NewAdvDif::new_assemble(arg_data_list &arg_data_v,const Nodedata *nodedata,
			     const Dofmap *dofmap,const char *jobinfo,
			     const ElementList &elemlist,
			     const TimeData *time_data) {

  int ierr=0;

  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_prof);

  int nelprops,nel,ndof;
  elem_params(nel,ndof,nelprops);
  int nen = nel*ndof;

  int locdof,kldof,lldof;
  char *value;

  SGETOPTDEF(int,double_layer,0);

  // allocate local vecs
  int kdof, kloc, node, jdim, ipg, nel2;
  FastMat2 veccontr(2,nel,ndof),xloc(2,nel,ndim),
    locstate(2,nel,ndof),vecc2;

  nen = nel*ndof;
  FastMat2 matloc(4,nel,ndof,nel,ndof),Jaco(2,ndimel,ndim),S(1,ndim),
    flux(1,ndof),load(1,ndof),hfilm(ndof);

  // Gauss Point data
  //o Type of element geometry to define Gauss Point data
  NGETOPTDEF_S(string,geometry,cartesian2d);
  GPdata gp_data(geometry.c_str(),ndim,nel,npg,GP_FASTMAT2);
  
  double detJ;
  FastMat2 state_in,state_out,u_in,u_out;

  double hfilm;
  if (double_layer) {
    PETSCFEM_ASSERT(nel % 2 ==0,"Number of nodes per element has to be even for "
		    "double_layer mode");
    nel2=nel/2;
    // state vector in both sides of the layer
    u_in.resize(2,nel2,ndof);
    u_out.resize(2,nel2,ndof);
    h_film_fun.init();
    ierr = get_double(thash,"hfilm",hfilm.Store(),0,ndof);
    CHKERRA(ierr);
    load=0;
    ierr = get_double(thash,"load",load.Store(),1,ndof); CHKERRA(ierr);
    locsta1.ReSize(nel2,ndof),locsta2.ReSize(nel2,ndof);
    xloc.ReSize(nel2,ndim);
    vecc2.ReSize(nel2,ndof);
  } else {
    ierr = get_double(thash,"load",load.Store(),0,ndof); CHKERRA(ierr);
    nel2=nel;
  }
  
  int ielh=-1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    ielh++;

    // Load local node coordinates in local vector
    for (kloc=0; kloc<nel2; kloc++) {
      node = ICONE(k,kloc);
      for (jdim=0; jdim<ndim; jdim++) {
	xloc(kloc+1,jdim+1) = NODEDATA(node-1,jdim);
      }
    }

    matloc = 0;
    veccontr = 0;
    locstate << &(LOCST(ielh,0,0));

    if (double_layer) {
      locsta1 = locstate.Rows(1,nel2);
      locsta2 = locstate.Rows(nel2+1,nel);
    }

#define DSHAPEXI (gp_data.dshapexi[ipg])
#define SHAPE    (gp_data.shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    // loop over Gauss points

    for (ipg=0; ipg<npg; ipg++) {

      Jaco = DSHAPEXI * xloc;

      if (ndimel==ndim) {
	detJ = mydet(Jaco);
      } else if (ndimel==ndim-1) {
	detJ = mydetsur(Jaco,S);
	S = -S; ; // fixme:= This is to compensate a bug in mydetsur
      } else {
	PFEMERRQ("Not considered yet ndimel<ndim-1\n");
      }

      if (!double_layer) {
	flux=load;
      } else {
	u1=SHAPE*locsta1;
	u2=SHAPE*locsta2;
	flux = (hfilm*(u2-u1)).t()+load;
      }
      
      double wpgdet = detJ*WPG;
      if (double_layer) {
	vecc2 = wpgdet * SHAPE.t() * flux;
	veccontr.Rows(1,nel2) += vecc2;
	veccontr.Rows(nel2+1,nel) -= vecc2;
      } else {
	veccontr  += wpgdet * SHAPE.t() * flux;
      }

    }

    if (ijob==COMP_VEC || ijob==COMP_FDJ || ijob==COMP_FDJ_PROF)
      veccontr >> &(RETVAL(ielh,0,0));
    if (ijob==COMP_MAT || ijob==COMP_MAT_PROF)
      matloc >> &(RETVALMAT(ielh,0,0,0,0));
    
  }
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
