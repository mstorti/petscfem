//__INSERT_LICENSE__
//$Id: metisprt.cpp,v 1.16 2002/05/06 14:07:33 mstorti Exp $

#include "fem.h"
#include "utils.h"
#include "util2.h"
#include "readmesh.h"
#include "idmap.h"
extern "C" {
#include <metis.h>
}

#include "elemset.h"
//#include "libretto.h"
#include <libretto/darray.h>
#include <libretto/autostr.h>
#include <libretto/autobuf.h>
#include <regex.h>
// STL components
#include <set>
#include <vector>
#include <map>
#include <cassert>
#include <queue>
#include "getprop.h"

#define ICONE(j,k) (icone[nel*(j)+(k)]) 

using namespace std;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "find_elem"
void find_elem(int elem,const int *nelemsetptr,int nelemsets,Mesh *mesh,
	       Elemset *& elemset,int &locel) {
  // return the elemset pointer and local elemset number of
  // a global elemset number
  for (int ielset=0; ielset<nelemsets; ielset++) {
    // Look for which elemset the element belongs
    if (elem<nelemsetptr[ielset+1]) {
      elemset  = *(Elemset **)da_ref(mesh->elemsetlist,ielset);
      locel = elem-nelemsetptr[ielset];
      return;
    }
  }
  assert(0);
}

// return a pointer to the connectivities of a particular element,
// given the global element number
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "elem_connectivities"
void elem_connectivities(int elem,Mesh *mesh,const int *nelemsetptr,
			 int nelemsets,
			 const int *&elem_icone,int &nel) {

  // iel:= the number of elemeset to which the element belongs
  // locel:= the element number local to the elemset
  int iel,locel,*icone;
  Elemset *elemset;
  
  find_elem(elem,nelemsetptr,nelemsets,mesh,elemset,locel);
  nel = elemset->nel;
  icone = elemset->icone;
  elem_icone = &ICONE(locel,0);
}

void  metis_part(int nelemfat,Mesh *mesh,
		 const int nelemsets,int *epart,
		 int *nelemsetptr,int *n2eptr,
		 int *node2elem,int size,const int myrank,
		 const int partflag,float *tpwgts,
		 int max_partgraph_vertices,
		 int iisd_subpart,int print_statistics) {

  Elemset *elemset;
  int *icone,nel,node,nvrtx,adjcount,j,elem,elemk,vrtx,
    visited,locel,k,jj,vrtxj,vrtxjj,p,ierr,nvsubdo;
  const int *elem_icone;
  double weight_scale=1.;
  
  // ngbrs:= Auxiliary queue that stores those nodes that are waiting
  // to be marked. 
  queue<int> ngbrs; 

  // Stores the graph `adjncy' in STL format
  vector< set<int> > adjncy_v;
  set<int>::iterator q,qe;
  set<int> subd_remap_aux;

  // nvrtx:= number of vertices used in graph partitioning
  // The factor 2 here avoids the case that only few elements are not
  // coalesced. In this case, for the last elements the computational
  // effort may be high, since the probability of getting a hole is
  // very low. 
  nvrtx = (nelemfat/2 > max_partgraph_vertices ?
	    max_partgraph_vertices : nelemfat);
  if (print_statistics && myrank==0) 
    printf("-- Graph partitioning statistics: ----------- \n"
	   "Using %d hraph vertices\n",nvrtx);

  // Create adjacency table for partitioning with Metis. In the
  // adjacency graph the nodes are elements or group of elments of the
  // FEM mesh. Two nodes of the graph (elements of the mesh) have are
  // linked if they share a node.

  // adjncy:= xadj:= graph described in CSR format (as defined in
  // Metis documentation)
  int *xadj = new int[nvrtx+1],*adjncy;
  // vwgt:= weights for the vertices (elements) of the graph
  int *vwgt = new int[nvrtx];
  // el2vrtx:= maps elements to vertices when coalescing
  int *el2vrtx = new int[nelemfat];
  // vpart:= partitioning of vertices
  int *vpart = new int[nvrtx];

  // Initialize 
  for (j=0; j<nelemfat; j++) el2vrtx[j] = -1;
  for (j=0; j<nvrtx; j++) vwgt[j] = 0;
  // Take `nvrtx' elements to coalesce the neighbor elements
  
  vrtx=0;
#if 0
  ngbrs.clear();
#else
  while (!ngbrs.empty()) ngbrs.pop();
#endif

  if (nelemfat!=nvrtx) {
    // Take nvrtx `seeds' for growing the groups 
    while(1) {
      elem = irand(nelemfat);
      if (el2vrtx[elem]!=-1) continue;
      el2vrtx[elem]=vrtx++;
      find_elem(elem,nelemsetptr,nelemsets,mesh,elemset,locel);
      vwgt[vrtx] += int(elemset->weight()/weight_scale);
      ngbrs.push(elem);
      if (vrtx==nvrtx) break;
    }
    visited = nvrtx;
  } else {
    // Take a group for each element
    for (elem=0; elem<nelemfat; elem++) {
      el2vrtx[elem] = elem;
      find_elem(elem,nelemsetptr,nelemsets,mesh,elemset,locel);
      vwgt[elem] += int(elemset->weight()/weight_scale);
    }      
    visited=nelemfat;
  }
  
  // Start coalescing neighbor elements until all elements are
  // assigned a vertex
  while (1) {
    if (ngbrs.empty()) {
      // either we have visited all the elements or there are
      // disconnected parts of the domain.

      // All elements have been visited
      if (visited==nelemfat) break;
      
      // Disconnected parts. Find next element not visited.
      for (k=0; k<nel; k++) {
	if (el2vrtx[k]<0) {
	  ngbrs.push(k);
	  // Put it in a random vertex
	  vrtx = irand(0,nvrtx)-1;
	  el2vrtx[k] = vrtx;
	  find_elem(k,nelemsetptr,nelemsets,mesh,elemset,locel);
	  vwgt[vrtx] += int(elemset->weight()/weight_scale);
	  break;
	}
      }
    }

    // Take the first element
    elem = ngbrs.front();
    ngbrs.pop();
    vrtx = el2vrtx[elem];
    assert(vrtx>=0);
    // Get its connectivities
    elem_connectivities(elem,mesh,nelemsetptr,nelemsets,elem_icone,nel);
    // Loop over all nodes of the element
    for (k=0; k<nel; k++) {
      node = elem_icone[k];
      // Loop over all elements connected to the node
      for (jj=n2eptr[node-1]; jj<n2eptr[node]; jj++) {
	elemk = node2elem[jj];
	// Already visited
	if (el2vrtx[elemk]>=0) continue;
	// Assign to group
	el2vrtx[elemk]=vrtx;
	visited++;
	// Find elemset and add weight to vertex weight (for load
	// balancing) 
	find_elem(elemk,nelemsetptr,nelemsets,mesh,elemset,locel);
	vwgt[vrtx] += int(elemset->weight()/weight_scale);
	// Put in queue
	ngbrs.push(elemk);
      }
    }
  }

#if 0
  if (myrank==0) {
    for (elem=0; elem < nelemfat; elem++) {
      printf("el2vrtx[%d] = %d\n",elem,el2vrtx[elem]);
    }
  }
#endif

  // mark:= auxiliary vector that flags if an element has been marked
  // already as a linked node in the graph
#if 0
  int *mark = new int[nelemfat];
  for (int ielgj=0; ielgj<nelemfat; ielgj++) mark[ielgj]=-1;
#endif
        
  adjncy_v.resize(nvrtx);

  // Define graph descriptors 
  // Initializa graph ptrs
  for (int ielset=0; ielset<nelemsets; ielset++) {
    elemset  = *(Elemset **)da_ref(mesh->elemsetlist,ielset);
    if (!elemset->isfat) continue;
    //if (myrank==0) printf("elemento: %d, conec: %d\n",ielset,elemset);
    icone = elemset->icone;
    nel = elemset->nel;
    const int *conn; int nell;
    for (int iel=0; iel<elemset->nelem; iel++) {
      int ielgj = nelemsetptr[ielset]+iel;
      vrtxj = el2vrtx[ielgj];
      nell = elemset->real_nodes(iel,conn);
      for (int iloc=0; iloc<nell; iloc++) {
	node = conn[iloc];
	
	// loop over all the elements connected to this node
	for (int jj=n2eptr[node-1]; jj<n2eptr[node]; jj++) {
	  int ielgjj = node2elem[jj];
	  if (ielgjj==ielgj) continue;
	  vrtxjj = el2vrtx[ielgjj];
	  adjncy_v[vrtxj].insert(vrtxjj);
	  // Just in case the user doesn't return a symmetric graph
	  adjncy_v[vrtxjj].insert(vrtxj);
	}
      }
    }
  }

  xadj[0] = 0;
  for (vrtxj=0; vrtxj<nvrtx; vrtxj++) 
    xadj[vrtxj+1] = xadj[vrtxj] + adjncy_v[vrtxj].size();

  // Once computed the adjncy size, it is created.
  // Computes statistics
  vector<int> vrtx_count;
  adjncy = new int[xadj[nvrtx]];
  for (vrtxj=0; vrtxj<nvrtx; vrtxj++) {
    set<int> &adj = adjncy_v[vrtxj];
    int e, adjs = adj.size();
    // Number of connected nodes may be zero (somewhat strange....)
    e = (adjs > 0 ? int(floor(log(adjs)/log(2.0)+1e-5)) : 0 );
    if (vrtx_count.size() <= e) vrtx_count.resize(e+1,0);
    vrtx_count[e]++;
    qe = adj.end();
    p = xadj[vrtxj];
    for (q=adj.begin(); q!=qe; q++) adjncy[p++] = *q;
  }

  if (print_statistics && myrank==0) {
    printf("Neighbor statistics for element graph:\n");
    int nne = 1;
    for (int e=0; e<vrtx_count.size(); e++) {
      printf("%5d graph vertices with %7d <= neighbors < %7d\n",
	     vrtx_count[e],nne,2*nne);
      nne *= 2;
    }
  }
  vrtx_count.clear();
      
#if 0
  // print the graph
  for (int ielgj=0; ielgj<nvrtx; ielgj++) {
    printf("%d: ",ielgj);
    for (int jj=xadj[ielgj]; jj<xadj[ielgj+1]; jj++) 
      printf(" %d",adjncy[jj]);
    printf(", weight: %d\n",vwgt[ielgj]);
  }
#endif

  int options=0,edgecut,numflag=0,wgtflag=2;
  // nvsubdo:= Nmber of `virtual' subdomains
  nvsubdo = size * iisd_subpart;
  float *tpwgts_d = new float[nvsubdo];
  // subd2proc:= this maps each subdomain to one processor. Initially
  // we put them random to avoid neighbor subdomains. The correct
  // solution should be to apply a colouring technique. Or better, to
  // repartition each subdomain.
  int *subd2proc = new int[nvsubdo];
  // subd_remap_aux:= An auxiliary STL set to obtain a random
  // permutation. First fill the set with `1,2,...,nvsubdo' and then
  // make a `random_pop' until it is void.
  for (j=0; j<nvsubdo; j++) subd_remap_aux.insert(j);
  for (j=0; j<nvsubdo; j++) subd2proc[j] = random_pop(subd_remap_aux);

  for (j=0; j<nvsubdo; j++) 
    tpwgts_d[j] = tpwgts[subd2proc[j]/iisd_subpart]/float(iisd_subpart);
  if (size*iisd_subpart > 1) {
    if (myrank==0) {
      if (partflag==0) {
	if (myrank==0) printf("METIS partition - partflag = %d\n",partflag);
	METIS_WPartGraphKway(&nvrtx,xadj,adjncy,vwgt, 
			     NULL,&wgtflag,&numflag,&nvsubdo, 
			     tpwgts_d,&options,&edgecut,vpart);
      
      } else { // partflag=2
	assert(0);
#if 0
	if (myrank==0) printf("Neighbor Partition - partflag = %d\n",partflag);
	int *mnnel = new int [size+1];
	for (int nnod=0; nnod<size; nnod++) {
	  mnnel[0]=0;
	  mnnel[nnod+1] = mnnel[nnod]+int(nelemfat*tpwgts[nnod]);
	  if (mnnel[size] < nelemfat) mnnel[size] = nelemfat;
	}
	// hago una busqueda de ngbrs por capas en adjncy y marco los
	// elementos (nodos del grafo) a los que le asigne una particion
	for (int imar=0; imar < nelemfat; imar++) {
	  vpart[imar]=0;
	}
	ngbrs.push(0);
	int counter=0;
	int inod=1;
	while (ngbrs.size()>0) {
	  int nods=ngbrs.front();
	  ngbrs.pop();
	  for (int vecnod_c=xadj[nods]; vecnod_c<xadj[nods+1];
	       vecnod_c++) {
	    int vecnod = adjncy[vecnod_c];
	    if (vecnod==nods) continue;
	    if (vpart[vecnod]!=0) continue;
	    ngbrs.push(vecnod);
	    vpart[vecnod]=inod; 
	  }
	  ++counter;
	  if (counter>mnnel[inod]) inod++;
	  vpart[nods]=inod;
	}
	for (int cnod=0; cnod<nelemfat; cnod++) {
	  vpart[cnod]=vpart[cnod]-1;
	}
#endif
      } // partflag==2
      
    } // if myrank==0
    // Remap subdomains to domains
    

    // Broadcast partitioning to other nodes
    ierr = MPI_Bcast(vpart,nvrtx,MPI_INT,0,PETSC_COMM_WORLD);
  } else {
    for (int jj=0; jj<nvrtx; jj++) 
      vpart[jj]=0;
  }
  for (elem=0; elem<nelemfat; elem++) 
    epart[elem] = subd2proc[vpart[el2vrtx[elem]]]/iisd_subpart;
  delete[] adjncy;
  delete[] xadj;
  delete[] el2vrtx;
  delete[] vwgt;
  delete[] vpart;
  delete[] tpwgts_d ;
  delete[] subd2proc;

#if 0
  for (vrtxj=0; vrtxj<nvrtx; vrtxj++)
    printf("vrtxj %d, proc %d\n",vrtxj,vpart[vrtxj]);
  for (elem=0; elem<nelemfat; elem++)
    printf("elem %d, proc %d\n",elem,epart[elem]);
#endif

  assert(ngbrs.empty());

}
