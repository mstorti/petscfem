//__INSERT_LICENSE__
//$Id: embgath.cpp,v 1.7 2002/08/07 15:26:33 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/linkgraph.h>

#include "embgath.h"
extern Mesh *GLOBAL_MESH;
extern int MY_RANK,SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
const int Quad2Hexa::faces[][8] = {
  0,1,2,3,4,5,6,7,
  1,5,6,2,0,4,7,3,
  0,4,5,1,3,7,6,2};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Surf2Vol::map_mask(const int *map_fc,int *vicorow) {
  int nel_surf, nel_vol, nf = nfaces(nel_surf,nel_vol);
  int match=0;
  const int *fc, *vol;
  for (int f=0; f<nf; f++) {
    face(f,fc,vol);
    match = 1;
    for (int l=0; l<nel_surf; l++) {
      if (map_fc[l] != fc[l]) {
	match=0;
	break;
      }
    }
    if (match) break;
  }
  assert(match);
  vector<int> vicorow_c(nel_vol);
  for (int j=0; j<nel_vol; j++) vicorow_c[j] = vicorow[vol[j]];
  for (int j=0; j<nel_vol; j++) vicorow[j] = vicorow_c[j];
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Quad2Hexa::face(int j,const int *&fc,const int *&vol_ret) { 
  int spin,rota,k,m;
  int spin_map[] = {0,3,2,1};
  spin = modulo(j,2,&m);
  rota = modulo(m,4,&k);
  for (int l=0; l<4; l++) {
    int ll = modulo(l+rota,4);
    vol[l] = faces[k][ll];
    vol[l+4] = faces[k][ll+4];
  }
  if (spin) {
    for (int l=0; l<4; l++) {
      vol_r[l] = vol[4+spin_map[l]];
      vol_r[4+l] = vol[spin_map[l]];
    }
    vol_ret = vol_r;
  } else vol_ret = vol;

  for (int k=0; k<4; k++) 
    this_face[k] = (use_exterior_normal() ? 
		    vol_ret[spin_map[k]] : vol_ret[k]);
  fc = this_face;
}

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
		"embedded_gatherer: Couldn't find volume_elemset = %s\n",
		volume_elemset.c_str());
    PetscFinalize();
    exit(0);
  }

  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  //o Number of Gauss points.
  TGETOPTNDEF(thash,int,npg,none);
  // ierr = get_int(thash,"npg",&npg); CHKERRA(ierr);
  TGETOPTNDEF(thash,int,ndim,none); //nd
  //o Use exterior or interior normal
  TGETOPTDEF(thash,int,use_exterior_normal,1);

  int ndimel=ndim-1;
  assert(geometry=="quad2hexa");
  Quad2Hexa gp_data(geometry.c_str(),ndim,nel,npg,
		    GP_FASTMAT2,use_exterior_normal);
  Surf2Vol *sv_gp_data = &gp_data;

  int nel_surf, nel_vol;
  surface_nodes(nel_surf,nel_vol);
  assert(nel_surf>0 && nel_surf<=nel);
  assert(nel_vol <= nel); //
  assert(nel_vol <= vol_elem->nel);
  
  // Mark nodos on the surface
  int nnod = GLOBAL_MESH->nodedata->nnod;
  // surface:= is surface[k]==0 then k is not on the surface
  // if != 0 then surface[k] is the number of surface node +1
  vector<int> surface(nnod,0);
  // maps surface numbering (0 to surf_nodes-1) to global (0 to nnod-1)
  vector<int> srf2glb;
  for (int e=0; e<nelem; e++) {
    int *icorow = icone + nel*e;
    for (int j=0; j<nel_surf; j++) surface[icorow[j]-1]=1;
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
    for (int j=0; j<nel_vol; j++) {
      int node = icorow[j]-1;
      int snode = surface[node]-1;
      if (snode>=0) graph.add(snode,e);
    }
  }

  // For each surface element look for the corresponding
  // volume element that shares a face
  vector<int> mask(nel_vol);
  for (int e=0; e<nelem; e++) {
    int *icorow = icone + nel*e;
    LinkGraphRow row;
    assert(nel_surf>0);
    // Take list for first node
    int node = icorow[0]-1;
    int sf_node = surface[node]-1;
    graph.set_ngbrs(sf_node,row);
    LinkGraphRow::iterator q;
    int found=0;
    int *vicorow;
    for (q=row.begin(); q!=row.end(); q++) {
      int ve = *q; // the volume element
      vicorow = vol_elem->icone + vol_elem->nel*ve;
      for (int j=0; j<nel_vol; j++) mask[j]=-1;
      found=0;
      for (int j=0; j<nel_surf; j++) {
	int sf_node = icorow[j];
	for (int k=0; k<nel_vol; k++) {
	  if (vicorow[k]==sf_node) {
	    mask[j] = k;
	    break;
	  }
	}
	if (mask[j]==-1) break;
	found++;
      }
      if (found==nel_surf) break;
    }
    for (int j=0; j<nel; j++) 
    if (found!=nel_surf) {
      PetscPrintf(PETSC_COMM_WORLD,
		  "embedded_gatherer: Can't find matching volume element"
		  " to surface element %d\n",e);
      PetscFinalize();
      exit(0);
    }
    for (int j=0; j<nel_vol; j++) icorow[j] = vicorow[j];
    // Volume element was found, find map and
    sv_gp_data->map_mask(mask.begin(),icorow);
  }

#if 1
  if (MY_RANK==0) {
    printf("Surface element connectivities: \n");
    for (int e=0; e<nelem; e++) {
      int *icorow = icone + nel*e;
      printf("surf.el. %d:",e+1);
      for (int j=0; j<nel_vol; j++) printf("%d ",icorow[j]);
      printf("\n");
    }
  }
#endif
  PetscFinalize();
  exit(0);
 
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
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
