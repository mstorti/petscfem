//__INSERT_LICENSE__
//$Id: embgath.cpp,v 1.3 2002/08/06 16:19:53 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/linkgraph.h>

#include "embgath.h"
extern Mesh *GLOBAL_MESH;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int embedded_gatherer::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(gather);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void embedded_gatherer::initialize() {
  int ierr;
  //o The name of the friend volume elemset
  TGETOPTDEF_S(thash,string,volume_elemset,);
  assert(volume_elemset.size()>0);
  int nelemsets = da_length(GLOBAL_MESH->elemsetlist);
  Elemset *vol_elem = NULL;
  for (int k=0; k<nelemsets; k++) {
    Elemset *e = *(Elemset **) da_ref(GLOBAL_MESH->elemsetlist,k);
    if (!strcmp(e->name(),volume_elemset.c_str())) {
      vol_elem = e;
      break;
    }
  }
  if (!vol_elem) {
    PetscPrintf(PETSC_COMM_WORLD,
		"embedded_gatherer: Couldn't find volume_elemset = %s\n"
		volume_elemset.c_str());
    PetscFinalize();
    exit(0);
  }

  // Construct graph for volume elemset
  LinkGraph graph;
  int nnod = GLOBAL_MESH->nodedata->nnod;
  graph.set_chunk_size(10000);
  graph.init(nnod);

  // Construct node to element array for the volume elemset
  for (int e=0; e<vol_elem->nelem; e++) {
    int *icorow = vol_elem->icone+vol_elem->nel*e;
    for (j=0; j<nel; j++) graph.add(icorow[j]-1,e);
  }
  
  for (int e=0; e<nelem; e++) {
    int *icorow = vol_elem->icone+vol_elem->nel*e;
  }
  graph.clear();
}
  
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int embedded_gatherer::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
		      Dofmap *dofmap,const char *jobinfo,int myrank,
		      int el_start,int el_last,int iter_mode,
		      const TimeData *time) {
  assert(0); // not defined yet
}


#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ
