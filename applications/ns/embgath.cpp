//__INSERT_LICENSE__
//$Id: embgath.cpp,v 1.5 2002/08/06 20:49:27 mstorti Exp $

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

  int nel_surf, nel_vol;
  surface_nodes(nel_vol,nel_vol);
  assert(nel_surf>0 && nel_surf<=nel);
  assert(nel_vol == vol_elem->nel);
  
  // Mark nodos on the surface
  int nnod = GLOBAL_MESH->nodedata->nnod;
  // surface:= is surface[k]==0 then k is not on the surface
  // if != 0 then surface[k] is the number of surface node +1
  vector<int> surface;
  // maps surface numbering (0 to surf_nodes-1) to global (0 to nnod-1)
  vector<int> srf2glb;
  for (int e=0; e<nelem; e++) {
    int *icorow = icone + nel*e;
    for (j=0; j<nel_surf; j++) surface[icorow[j]-1]=1;
  }
  // Count surface nodes
  int surf_nodes = 0;
  for (int k=0; k<nnod; k++) {
    if (surface[k]) {
      surface[k] = ++surf_nodes;
      srf2glb.push_back(surf_nodes);
    }
  }

  // Construct graph for volume elemset
  LinkGraph graph;
  graph.set_chunk_size(10000);
  graph.init(surf_nodes);

  // Construct node to element array for the volume elemset
  for (int e=0; e<vol_elem->nelem; e++) {
    int *icorow = vol_elem->icone + vol_elem->nel*e;
    for (j=0; j<nel_vol; j++) {
      int node = icorow[j]-1;
      int snode = surface[node];
      if (snode) graph.add(snode,e);
    }
  }

  // This should change with geometry
//    int nfaces = 6; // Number of faces per volume element
//    int faces[][] = {
//      { 
//    };
  
  // For each surface element look for the corresponding
  // volume element that shares a face
  vector<int> mask(nel_vol);
  for (int e=0; e<nelem; e++) {
    int *icorow = icone + nel*e;
    LinkGraphRow &row;
    graph.set_ngbrs(sf_node,row);
    LinkGraphRow::iterator q;
    for (q=row.begin(); q!=row.end(); q++) {
      int ve = *q; // the volume element
      int *vicorow = vol_elem->icone + vol_elem->nel*ve;
      for (int j=0; j<nel_vol; j++) mask[j]=0;
      int found=1;
      for (int j=0; j<nel_surf; j++) {
	int sf_node = icorow[j];
	found=0;
	for (int k=0; k<nel_vol; k++) {
	  if (vicorow[k]==sf_node) {
	    mask[k] = j;
	    found = 1;
	    break;
	  }
	}
	if (found) break;
      }
      if (!found) {
	PetscPrintf(PETSC_COMM_WORLD,
		    "embedded_gatherer: Can't find matching volume element"
		    " to surface element %d\n",e);
	PetscFinalize();
	exit(0);
      }
 
    }
  }

  graph.clear();
  surface.clear();
  srf2glb.clear();
}
  
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int embedded_gatherer::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
		      Dofmap *dofmap,const char *jobinfo,int myrank,
		      int el_start,int el_last,int iter_mode,
		      const TimeData *time) {
  assert(0); // not defined yet
}

void visc_force_integrator::surface_nodes(int &nel_surf,int &nel_vol)=0 {
}

void visc_force_integrator
::set_pg_values(vector<double> &pg_values,FastMat2 &u,
		FastMat2 &uold,FastMat2 &grad_u, FastMat2 &grad_uold, 
		FastMat2 &xpg,FastMat2 &n,
		double wpgdet,double time) {
  assert(0); // not defined yet
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void visc_force_integrator
::surface_nodes(int &nel_surf,int &nel_vol) { 
  nel_surf=4; 
  nel_vol=8; 
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ
