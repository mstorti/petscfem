//__INSERT_LICENSE__
//$Id: gatherer.cpp,v 1.1 2002/03/16 22:24:17 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "gatherer.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int ns_volume_element::ask(char *,int &)"
int gatherer::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(gather);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "nsi_tet_les_fm2::assemble"
int gatherer::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
		      Dofmap *dofmap,const char *jobinfo,int myrank,
		      int el_start,int el_last,int iter_mode,
		      const TimeData *time) {

  GET_JOBINFO_FLAG(gather);
  assert(gather);

  //o Number of Gauss points.
  TGETOPTNDEF(thash,int,npg,none);
  //o Dimension od the embedding space
  TGETOPTNDEF(thash,int,ndim,none);
  //o Dimenson od the element
  TGETOPTNDEF(thash,int,ndimel,ndim-1); 
  //o Defines the geomtry of the element
  NGETOPTDEF_S(string,geometry,cartesian2d);

  GPdata gp_data(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);

#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)

  ja = 0;
  locst = arg_data_v[ja++].locst;
  locst2 = arg_data_v[ja++].locst;
  int options = arg_data_v[ja].options;
  vector<double> *values = arg_data_v[ja++].vector_assoc;
  int nvalues = values.size();
  vector<double> pg_values(nvalues);

  xloc.resize(2,nel,ndim);

  FastMat2 Jaco(ndimel,ndim),staten(2,nel,ndof), 
    stateo(2,nel,ndof),u_old(1.ndof),u(1,ndof);

  Time * time_c = (Time *)time;
  double t = time_c.time();
  
  // Initialize the call back functions
  init();

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  int ielh=-1,kloc,nel,ndof;
  for (int k=el_start; k<=el_last; k++) {

    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(k,kloc);
      xloc.ir(1,kloc+1).set(&NODEDATA(node-1,0));
    }

    staten.set(&(LOCST(ielh,0,0)));
    stateo.set(&(LOCST2(ielh,0,0)));

    // Let user do some things when starting with an element
    element_hook(k);

    for (ipg=0; ipg<npg; ipg++) {
      // Gauss point coordinates
      xpg.prod(SHAPE,xloc,-1,-1,1);
      // Jacobian master coordinates -> real coordinates
      Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);
     
      detJaco = Jaco.det();
      if (detJaco <= 0.) {
	printf("Jacobian of element %d is negative or null\n"
	       " Jacobian: %f\n",k,detJaco);
	PetscFinalize();
	exit(0);
      }
      wpgdet = detJaco*WPG;

      if (ndim==ndimel) {
	iJaco.inv(Jaco);
	detJaco = Jaco.det();
      } else if (ndimel==1) {
	// This allows to solve problems on streams like rivers or
	// ducts or advective problems on plane surfaces (not
	// implemented yet). Also, it could be used also for advective
	// problems on arbitrary surfaces (ndim=3 and ndimel=2) but I
	// don't know how to do that yet. (tensorial calculus...)
	detJaco = Jaco.norm_p_all(2);
	iJaco.setel(1./detJaco,1,1);
      }
      
      // Values of variables at Gauss point
      u.prod(SHAPE,u,-1,-1,1);
      u_old.set(SHAPE,stateo,-1,-1,1);

      set_pg_values(pg_values,u,uold,xpg,Jaco,wpgdet,t);
      if (options & VECTOR_ADD) {
	for (int j=0; j<nvalues; j++) {
	  values[j] += wpgdet*pg_values[j];
	}
      } else assert(0);
    }
  }  
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ
