//__INSERT_LICENSE__
//$Id: gatherer.cpp,v 1.4 2005/02/21 18:50:28 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "./gatherer.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int gatherer::ask"
int gatherer::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(gather);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "gatherer::assemble"
int gatherer::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
		      Dofmap *dofmap,const char *jobinfo,int myrank,
		      int el_start,int el_last,int iter_mode,
		      const TimeData *time) {

  int ierr;

  GET_JOBINFO_FLAG(gather);
  assert(gather);

  //o Position in gather vector
  TGETOPTDEF(thash,int,gather_pos,0);
  //o How many gather values will be contributed by this elemset
  TGETOPTDEF_ND(thash,int,gather_length,0);
  //o Number of Gauss points.
  TGETOPTNDEF(thash,int,npg,none);
  //o Dimension of the embedding space
  TGETOPTNDEF(thash,int,ndim,none);
  //o Dimenson of the element
  TGETOPTDEF(thash,int,ndimel,ndim-1); 
  //o Defines the geomtry of the element
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);

  GPdata gp_data(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);

#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)

  // get number of fields per node (constant+variables)
  int nu=nodedata->nu;

  int ja = 0;
  double *locst = arg_data_v[ja++].locst;
  double *locst2 = arg_data_v[ja++].locst;
  int options = arg_data_v[ja].options;
  vector<double> *values = arg_data_v[ja++].vector_assoc;
  int nvalues = values->size();
  assert(gather_pos+gather_length <= nvalues); // check that we don't put values
				   // beyond the end of global vector
				   // `values'
  vector<double> pg_values(gather_length);

  FastMat2 xloc(2,nel,ndim);

  FastMat2 Jaco(2,ndimel,ndim),staten(2,nel,ndof), 
    stateo(2,nel,ndof),u_old(1,ndof),u(1,ndof),
    n(1,ndim),xpg(1,ndim);

  Time * time_c = (Time *)time;
  double t = time_c->time();
  
  // Initialize the call back functions
  init();

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  int ielh=-1,kloc;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    FastMat2::reset_cache();
    ielh++;

    for (kloc=0; kloc<nel; kloc++) {
      int node = ICONE(k,kloc);
      xloc.ir(1,kloc+1).set(&NODEDATA(node-1,0));
    }
    xloc.rs();

    staten.set(&(LOCST(ielh,0,0)));
    stateo.set(&(LOCST2(ielh,0,0)));

    // Let user do some things when starting with an element
    element_hook(k);

    for (int ipg=0; ipg<npg; ipg++) {
      // Gauss point coordinates
      xpg.prod(SHAPE,xloc,-1,-1,1);
      // Jacobian master coordinates -> real coordinates
      Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);

      double detJaco;
      if (ndimel==ndim) {
	detJaco = Jaco.det();
      } else if (ndimel==ndim-1) {
	detJaco = Jaco.detsur(&n);
	n.scale(1./detJaco);
	n.scale(-1.);		// fixme:= This is to compensate a bug in mydetsur
      }
      if (detJaco <= 0.) {
	detj_error(detJaco,k);
	set_error(1);
      }
      double wpgdet = detJaco*WPG;

      // Values of variables at Gauss point
      u.prod(SHAPE,staten,-1,-1,1);
      u_old.prod(SHAPE,stateo,-1,-1,1);

      set_pg_values(pg_values,u,u_old,xpg,n,wpgdet,t);
      if (options & VECTOR_ADD) {
	for (int j=0; j<gather_length; j++) {
	  (*values)[gather_pos+j] += pg_values[j];
	}
      } else assert(0);
    }
  }  
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  return 0;
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ
