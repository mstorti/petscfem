//__INSERT_LICENSE__
//$Id: metisprt.cpp,v 1.2 2001/11/11 12:53:42 mstorti Exp $
 
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

void  metis_part(int nelemfat,Mesh *mesh,
		 const int nelemsets,int *vpart,
		 int *nelemsetptr,int *n2eptr,
		 int *node2elem,int size,const int myrank,
		 const int partflag,float *tpwgts,
		 int max_partgraph_vertices) {

  Elemset *elemset;
  int *icone,nel,node,nvertx;
  queue<int> vecinos; 

  // nvertx:= number of vertices used in graph partitioning
  nvertx = (nelemfat>max_partgraph_vertices ?
	    max_partgraph_vertices : nelemfat);

  // Create adjacency table for partitioning with Metis. In the adjacency
  // graph the nodes are the elements of the FEM mesh. Two nodes of
  // the graph (elements of the mesh) have are linked if they share a
  // node. 

  // adjncy:= xadj:= graph desdcribed in CSR format (as defined in
  // Metis documentation)
  int *xadj = new int[nvertx+1];

  // mark:= auxiliary vector that flags if an element has been marked
  // already as a linked node in the graph
  int *mark = new int[nelemfat];
  for (int ielgj=0; ielgj<nelemfat; ielgj++) 
    mark[ielgj]=-1;
    
    // Define graph descriptors 
  xadj[0]=0;
  for (int ielset=0; ielset<nelemsets; ielset++) {
    elemset  = *(Elemset **)da_ref(mesh->elemsetlist,ielset);
    //if (myrank==0) printf("elemento: %d, conec: %d\n",ielset,elemset);
    icone = elemset->icone;
    nel = elemset->nel;
    for (int iel=0; iel<elemset->nelem; iel++) {
      int ielgj = nelemsetptr[ielset]+iel;
      xadj[ielgj+1]=xadj[ielgj];
      // loop over the connected nodes
      for (int iloc=0; iloc<elemset->nel; iloc++) {
	node = ICONE(iel,iloc);
	
	// loop over all the elements connected to this node
	for (int jj=n2eptr[node-1];
	     jj<n2eptr[node]; jj++) {
	  int ielgjj = node2elem[jj];
	  if (ielgjj==ielgj) continue;

	  // check if the element has been already loaded
	  if (mark[ielgjj]!=ielgj) {
	    mark[ielgjj]=ielgj;
	    xadj[ielgj+1]++;
	    // printf("adding %d %d\n",ielgj,ielgjj);
	  }
	}
      }
    }
  }
  // Once computed the adjncy size, it is created.
  int adycont=0;
  int *adjncy = new int[xadj[nelemfat]]; // fixme:= identifier adjncy is overloaded
  for (int ielgj=0; ielgj<nelemfat; ielgj++) 
    mark[ielgj]=-1;
  // Define graph descriptors 
  xadj[0]=0;
  for (int ielset=0; ielset<nelemsets; ielset++) {
    elemset  = *(Elemset **)da_ref(mesh->elemsetlist,ielset);
    icone = elemset->icone;
    nel = elemset->nel;
    for (int iel=0; iel<elemset->nelem; iel++) {
      int ielgj = nelemsetptr[ielset]+iel;
      xadj[ielgj+1]=xadj[ielgj];
      // loop over the connected nodes
      for (int iloc=0; iloc<elemset->nel; iloc++) {
	node = ICONE(iel,iloc);
	
	// loop over all the elements connected to this node
	for (int jj=n2eptr[node-1];
	     jj<n2eptr[node]; jj++) {
	  int ielgjj = node2elem[jj];
	  if (ielgjj==ielgj) continue;

	  // check if the element has been already loaded
	  if (mark[ielgjj]!=ielgj) {
	    adjncy[adycont]=ielgjj;
	    mark[ielgjj]=ielgj;
	    xadj[ielgj+1]++;
	    adycont++;
	  }
	}
      }
    }
  }
  delete[] mark;

#if 0
  // print the graph
  for (int ielgj=0; ielgj<nelemfat; ielgj++) {
    printf("%d: ",ielgj);
    for (int jj=xadj[ielgj]; jj<xadj[ielgj+1]; jj++) 
      printf(" %d",adjncy[jj]);
    printf("\n");
  }
#endif

  int options=0,edgecut,numflag=0,wgtflag=0;
  if (size>1) {
    if (partflag==0) {
      if (myrank==0) printf("METIS partition - partflag = %d\n",partflag);
      METIS_WPartGraphKway(&nelemfat,xadj,adjncy,NULL, 
			   NULL,&wgtflag,&numflag,&size, 
			   tpwgts,&options,&edgecut,vpart);
    } else { // partflag=2
      if (myrank==0) printf("Neighbor Partition - partflag = %d\n",partflag);
      int *mnnel = new int [size+1];
      for (int nnod=0; nnod<size; nnod++) {
	mnnel[0]=0;
	mnnel[nnod+1] = mnnel[nnod]+int(nelemfat*tpwgts[nnod]);
	if (mnnel[size] < nelemfat) mnnel[size] = nelemfat;
      }
      // hago una busqueda de vecinos por capas en adjncy y marco los
      // elementos (nodos del grafo) a los que le asigne una particion
      for (int imar=0; imar < nelemfat; imar++) {
	vpart[imar]=0;
      }
      vecinos.push(0);
      int counter=0;
      int inod=1;
      while (vecinos.size()>0) {
	int nods=vecinos.front();
	vecinos.pop();
	for (int vecnod_c=xadj[nods]; vecnod_c<xadj[nods+1];
	     vecnod_c++) {
	  int vecnod = adjncy[vecnod_c];
	  if (vecnod==nods) continue;
	  if (vpart[vecnod]!=0) continue;
	  vecinos.push(vecnod);
	  vpart[vecnod]=inod; 
	}
	++counter;
	if (counter>mnnel[inod]) inod++;
	vpart[nods]=inod;
      }
      for (int cnod=0; cnod<nelemfat; cnod++) {
	vpart[cnod]=vpart[cnod]-1;
      }
    }
    // length(vpart) = number of nelemfat in the graph
  } else {
    for (int jj=0; jj<nelemfat; jj++) 
      vpart[jj]=0;
  }
  delete[] adjncy;
  delete[] xadj;
  assert(vecinos.empty());

}
