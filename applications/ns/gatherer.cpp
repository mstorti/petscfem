//__INSERT_LICENSE__
//$Id: gatherer.cpp,v 1.2 2002/03/17 03:53:40 mstorti Exp $

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

  int ierr;

  GET_JOBINFO_FLAG(gather);
  assert(gather);

  //o Number of Gauss points.
  TGETOPTNDEF(thash,int,npg,none);
  //o Dimension od the embedding space
  TGETOPTNDEF(thash,int,ndim,none);
  //o Dimenson od the element
  TGETOPTNDEF(thash,int,ndimel,ndim-1); 
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
  vector<double> pg_values(nvalues);

  FastMat2 xloc(2,nel,ndim);

  FastMat2 Jaco(ndimel,ndim),staten(2,nel,ndof), 
    stateo(2,nel,ndof),u_old(1,ndof),u(1,ndof);

  Time * time_c = (Time *)time;
  double t = time_c->time();
  
  // Initialize the call back functions
  init();

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  int ielh=-1,kloc,nel,ndof;
  for (int k=el_start; k<=el_last; k++) {

    for (kloc=0; kloc<nel; kloc++) {
      int node = ICONE(k,kloc);
      xloc.ir(1,kloc+1).set(&NODEDATA(node-1,0));
    }

    staten.set(&(LOCST(ielh,0,0)));
    stateo.set(&(LOCST2(ielh,0,0)));

    // Let user do some things when starting with an element
    element_hook(k);

    for (int ipg=0; ipg<npg; ipg++) {
      // Gauss point coordinates
      xpg.prod(SHAPE,xloc,-1,-1,1);
      // Jacobian master coordinates -> real coordinates
      Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);
     
      double detJaco = Jaco.det();
      if (detJaco <= 0.) {
	printf("Jacobian of element %d is negative or null\n"
	       " Jacobian: %f\n",k,detJaco);
	PetscFinalize();
	exit(0);
      }
      double wpgdet = detJaco*WPG;

      if (ndim==ndimel) {
	detJaco = Jaco.det();
      } else if (ndimel==1) {
	// This allows to solve problems on streams like rivers or
	// ducts or advective problems on plane surfaces (not
	// implemented yet). Also, it could be used also for advective
	// problems on arbitrary surfaces (ndim=3 and ndimel=2) but I
	// don't know how to do that yet. (tensorial calculus...)
	detJaco = Jaco.norm_p_all(2);
      }
      
      // Values of variables at Gauss point
      u.prod(SHAPE,u,-1,-1,1);
      u_old.prod(SHAPE,stateo,-1,-1,1);

      set_pg_values(pg_values,u,u_old,xpg,Jaco,wpgdet,t);
      if (options & VECTOR_ADD) {
	for (int j=0; j<nvalues; j++) {
	  (*values)[j] += wpgdet*pg_values[j];
	}
      } else assert(0);
    }
  }  
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void force_integrator::set_pg_values(vector<double> &pg_values,FastMat2 &u,
				     FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
				     double wpgdet,double time) {
  assert(Jaco.dim(1)==3);
  assert(Jaco.dim(2)==2);
  
  double v;
  v = Jaco.get(1,2)*Jaco.get(2,3) - Jaco.get(2,2)*Jaco.get(1,3);
  normal.setel(v,1);
  v = Jaco.get(1,3)*Jaco.get(2,1) - Jaco.get(2,3)*Jaco.get(1,1);
  normal.setel(v,2 );
  v = Jaco.get(1,1)*Jaco.get(2,2) - Jaco.get(2,1)*Jaco.get(1,2);
  normal.setel(v,3);

  v = sqrt(normal.sum_square_all());
  normal.scale(1./v);

  pg_values[0] += wpgdet*u.get(4)*normal.get(1);
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ
