//__INSERT_LICENSE__
//$Id: measperf.cpp,v 1.7 2004/07/28 22:06:18 mstorti Exp $

#include "fem.h"
#include <set>
#include "utils.h"
#include "getprop.h"
#include "elemset.h"
#include "idmap.h"
#include "dofmap.h"
#include "arglist.h"

// iteration modes
#define NOT_INCLUDE_GHOST_ELEMS 0
#define INCLUDE_GHOST_ELEMS 1
extern int MY_RANK,SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int measure_performance_fun()"
int measure_performance_fun(Mesh *mesh,arg_list argl,
			Dofmap *dofmap,const char *jobinfo,const TimeData
			*time) {
  PetscPrintf(PETSC_COMM_WORLD,
	      "----\n Measuring performance for jobinfo \"%s\"...\n",jobinfo);
  if (SIZE>1) 
    PetscPrintf(PETSC_COMM_WORLD,
		"Warning: measuring performance with size>1 !!!\n");

  Chrono chrono;

  double wnelem=0.;
  //o Number of times the loop is executed when measuring
  // performance. 
  int ierr;
  TGETOPTNDEF(mesh->global_options,int,measure_performance_loop_length,10);
#define NLOOP measure_performance_loop_length
  
  for (int ielset=0; ielset<da_length(mesh->elemsetlist); ielset++) {
    Elemset *elemset  = *(Elemset **)da_ref(mesh->elemsetlist,ielset);
    
    //o Weight for computing performance
    TGETOPTDEF(elemset->thash,double,measure_performance_weight,1.);
    TGETOPTDEF_S(elemset->thash,string,name,""); //nd 

    wnelem += elemset->nelem * measure_performance_weight;
    PetscPrintf(PETSC_COMM_WORLD,
		" Elemset type: %s, name: %s,\n    ptr: %p, "
		"elem: %d, weight per elem.: %.4g\n",
		elemset->type, (name == "" ? "<anonymous>" : name.c_str()),
		elemset, elemset->nelem,
		measure_performance_weight);
  }
  PetscPrintf(PETSC_COMM_WORLD,
	      "Total (weighted) number of elements: %.4f\n",
	      wnelem);

  double tsum=0., t2sum=0.,tmax, tmin,mean,max,min,dev,t;
  for (int jjj=0; jjj<NLOOP; jjj++) {
    chrono.start();
    ierr = assemble(mesh,argl,dofmap,jobinfo,time); CHKERRA(ierr);
    t = chrono.elapsed();
    printf("[loop iter %d] elapsed %.4g\n",jjj+1,t);
    tsum += t;
    t2sum += t*t;
    if (jjj==0) {
      tmin=t;
      tmax=t;
    }
    if (t>tmax) tmax=t;
    if (t<tmin) tmin=t;
  }

  max = tmax/wnelem*1000.;
  min = tmin/wnelem*1000;
  mean = tsum/NLOOP/wnelem*1000;
  dev = sqrt(t2sum/NLOOP-tsum*tsum/NLOOP/NLOOP)/wnelem*1000;
  PetscPrintf(PETSC_COMM_WORLD,"Total %.4g, ntimes %d, nelems %.4f\n",
	 tsum, NLOOP, wnelem);
  PetscPrintf(PETSC_COMM_WORLD,
	      "Rate [sec/Kelems]: mean %.4g, min %.4g, "
	      "max %.4g, std.dev. %.4g\n",mean,min,max,dev);
  return 0;
}
#undef NLOOP

