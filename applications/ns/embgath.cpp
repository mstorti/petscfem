//__INSERT_LICENSE__
//$Id: embgath.cpp,v 1.2 2002/08/06 15:31:20 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

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
  int nelemsets = da_length(GLOBAL_MESH->elemsetlist);
  Elemset *vol_elem = NULL;
  for (int k=0; k<nelemsets; k++) {
    Elemset *e = *(Elemset **) da_ref(GLOBAL_MESH->elemsetlist,k);
    if (!strcmp(e->name(),volume_elemset.c_str())) {
      vol_elem = e;
      printf("I found it!!!\n");
      break;
    }
  }
  if (!vol_elem) printf("Couldn't find it!!!\n");
  PetscFinalize();
  exit(0);
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
