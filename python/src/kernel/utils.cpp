//$Id$

// STL components
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <cassert>
using namespace std;


#include "petscfem4py.h"
#include "Error.h"


typedef PF4PY_NAMESPACE::Error Error;


// PETSc-FEM components
#include <src/fem.h>
#include <src/elemset.h>
#include <src/dofmap.h>

extern void
metis_part(MPI_Comm comm, int nelemfat,Mesh *mesh,
	   const int nelemsets,int *vpart,
	   int *nelemsetptr,int *n2eptr,
	   int *node2elem,int size,const int rank,
	   const int partflag,float *tpwgts,
	   int max_partgraph_vertices,
	   int iisd_subpart,
	   int print_partitioning_statistics);

#undef  ICONE
#define ICONE(j,k) (icone[nel*(j)+(k)]) 

#undef  __FUNC__
#define __FUNC__ "mesh_setup"
static int 
mesh_setup(MPI_Comm comm, Mesh* mesh,
	   std::vector<int>&   nodepart,
	   std::vector<float>& weights)
{
  PetscFunctionBegin;
  
  assert(mesh != NULL);

  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  
  int nnod = mesh->nodedata->nnod;
  
  float* tpwgts    = &weights[0];
  int*   npart     = &nodepart[0];

  int ierr;
  int node, jdof,kdof;
  int nelem, nel, iele;
  int *icone;
  row_t row,col; row_t::iterator kndx;

  std::size_t nelemsets = da_length(mesh->elemsetlist);
  int nelemfat=0, is_any_fat=0;
  int *n2eptr, *nelemsetptr;


  for (std::size_t ielset=0; ielset<nelemsets; ielset++) {

    Elemset** eptr = (Elemset**) da_ref(mesh->elemsetlist, ielset);
    assert (eptr != NULL);
    Elemset* elemset  = *(Elemset **)da_ref(mesh->elemsetlist,ielset);
    assert (elemset != NULL);
    
    elemset->e1          = 0;
    elemset->e2          = 0;
    elemset->isfat       = 1;
    elemset->nelem_here  = 0;
    memset(elemset->elem_conne, 0, sizeof(int)*elemset->nel);
    memset(elemset->epart,      0, sizeof(int)*elemset->nelem);
    memset(elemset->epart2,     0, sizeof(int)*elemset->nelem);
    da_resize(elemset->ghost_elems, 0);
  }

  for (int k=0; k<nnod; npart[k++]=1);

  std::vector<int> _v_n2eptr(nnod+1, 0);
  std::vector<int> _v_nelemsetptr(nelemsets+1, 0);
  n2eptr = &_v_n2eptr[0];
  nelemsetptr = &_v_nelemsetptr[0]; // @@@

  // compute total number of elements in fat elemsets, nelemfat
  // also n2eptr[node] gets the total number of elements connected to
  // node and n2eptr[0]=0

  // numfat:= number of fat elemsets
  for (std::size_t ielset=0; ielset<nelemsets; ielset++) {
    Elemset * elemset  = *(Elemset **)da_ref(mesh->elemsetlist,ielset);
    int is_fat=1;
    ierr = get_int(elemset->thash,"is_fat",&is_fat,1); CHKERRQ(ierr);
    elemset->isfat=is_fat;
    nelemsetptr[ielset+1]=nelemsetptr[ielset];
    if (!is_fat) continue;
    nelem = elemset->nelem;
    icone = elemset->icone;
    nel = elemset->nel;
    is_any_fat=1;
    nelemfat += nelem;
    nelemsetptr[ielset+1]=nelemsetptr[ielset]+nelem;
    const int *conn; int nell;
    for (int iel=0; iel<nelem; iel++) {
      nell = elemset->real_nodes(iel,conn);
      for (int iloc=0; iloc<nell; iloc++) {
	// printf("adding to node %d\n",CONN(iel,iloc));
	n2eptr[conn[iloc]]++;
      }
    }
  }

  // raise error if not defined at least one fat elemset
  if (size > 1 && !is_any_fat)
    throw Error("Mesh: at least one elemset has to be set as fat "
		"(set option 'is_fat' to '1') "
		"if number of processors is greater than one");

  // cumulated sum 
  for (node=0; node<nnod; node++) 
    n2eptr[node+1] += n2eptr[node];
  // PetscPrintf(comm,"n2eptr[nnod]: %d\n",n2eptr[nnod]);

  // n2esize:= total number of entries in the node2elem table
  int n2esize = n2eptr[nnod];

  // initialize node2elem
  std::vector<int> _v_node2elem(n2esize, -1);
  int *node2elem = &_v_node2elem[0];

  // define node2elem pointer
  for (std::size_t ielset=0; ielset<nelemsets; ielset++) {
    Elemset* elemset  = *(Elemset **)da_ref(mesh->elemsetlist,ielset);
    if (!elemset->isfat) continue;
    nelem = elemset->nelem;
    // icone = elemset->icone;
    nel = elemset->nel;
    for (int iel=0; iel<nelem; iel++) {
      // iel:= global element number cumulated through the elemsetlist,
      // (only on fat elemsets)
      int ielg = nelemsetptr[ielset]+iel;
      const int *conn; int nell;
      nell = elemset->real_nodes(iel,conn);
      for (int iloc=0; iloc<nell; iloc++) {
	// search for the next free position on the array
	node=conn[iloc];
	int jpos;
	for (jpos=n2eptr[node-1]; jpos<n2eptr[node]; jpos++) 
	  if (node2elem[jpos]==-1) break;
	if (jpos==n2eptr[node]) {
	  printf("node %d, elems %d, range in node2elem %d - %d\n",
		 node,n2eptr[node]-n2eptr[node-1],n2eptr[node-1],
		 n2eptr[node]-1);
	  for (jpos=n2eptr[node-1]; jpos < n2eptr[node]; jpos++)
	    printf("%d ",node2elem[jpos]);
	  printf("\n");
	  assert(0);
	}
	node2elem[jpos]=ielg;
      }
    }
  }

  // partflag = 0 -> METIS partition;
  // partflag = 1 -> Hitchhiking partition;
  // partflag = 2 -> Neighbor partition. 
  // partflag = 3 -> random partition. 
  // partflag = 4 -> natural partition. 
  int partflag=0;
  //o Set partitioning method. May be set to  #metis# ,
  //  #hitchhiking# ,  #nearest_neighbor#  or  #random# .
  TGETOPTDEF_S(mesh->global_options,string,partitioning_method,metis);
  //o Print graph statistics
  TGETOPTDEF(mesh->global_options,int,print_partitioning_statistics,0);
#define INF INT_MAX
  //o Maximum number of vertices admissible while computing the
  // partitioning graph.
  TGETOPTDEF(mesh->global_options,int,max_partgraph_vertices,INF);
#undef INF
  // o Number of subpartitions inside each processor. 
  TGETOPTDEF(mesh->global_options,int,iisd_subpart,1); //nd
  iisd_subpart = 1; // In order to deactivate subpartitioning at the high level

  if (partitioning_method == string("metis")) {
    partflag = 0;
  } else if (partitioning_method == string("hitchhiking")) {
    partflag = 1;
  } else if (partitioning_method == string("nearest_neighbor")) {
    partflag = 2;
  } else if (partitioning_method == string("random")) {
    partflag = 3;
  } else if (partitioning_method == string("natural")) {
    partflag = 4;
  } else {
    throw Error("Mesh: partitioning method not known");
  }

  ////PetscPrintf(comm,"Starts partitioning.\n"); 
  //std::vector<int> _v_vpart(nelemfat, 0);
  //int *vpart = &_v_vpart[0];
  int *vpart = new int[nelemfat];
  memset(vpart, 0, sizeof(int)*nelemfat);

  if (partflag==0 || partflag==2) {
    
    metis_part(comm,nelemfat,mesh,nelemsets,vpart,
	       nelemsetptr,n2eptr,node2elem,size,rank,
	       partflag,tpwgts,max_partgraph_vertices,
	       iisd_subpart,print_partitioning_statistics);

  } else if (partflag==1) {
    std::vector<int> _v_mnel(size+1, 0);
    int *mnel = &_v_mnel[0];
    if (size>1) {
      ///if (rank==0) printf("Hitchicking partition, partflag = %d\n",partflag);
      for (int nproc=0; nproc < size; nproc++) {
	mnel[0] = 0;
	mnel[nproc+1] = mnel[nproc]+int(nelemfat*tpwgts[nproc]);
	if (mnel[size] < nelemfat) mnel[size] = nelemfat-1;
	for (int jelem=mnel[nproc]; jelem < mnel[nproc+1]+1; jelem++) { 
	  vpart[jelem]=nproc;
	}
      }
    } else {
      for (int jjy=0; jjy<nelemfat; jjy++) 
	vpart[jjy]=0;
    }	

  } else if (partflag==3 || partflag==4) {
    // random partitioning 
    if (rank==0) {
      for (int j=0; j < nelemfat; j++) {
	vpart[j] = int(drand()*double(size));
	if (partflag==4)
	  vpart[j] = int(double(j)/double(nelemfat)*double(size));
	// Just in case random functions are too close to the limits
	// (In theory, rand() should not touch the limits 0, RAND_MAX)
	if (vpart[j]>=size) vpart[j]=size-1;
	if (vpart[j]<0) vpart[j]=0;
	// printf("element %d in [%d]\n",j+1,vpart[j]);
      }
    }
    ierr = MPI_Bcast (vpart,nelemfat,MPI_INT,0,comm); CHKERRQ(ierr);
  } else assert(0); // something went wrong

  // nelem_part:= nelem_part[proc] is the number of elements in
  // processor proc. 
  std::vector<int> _v_nelem_part(size, 0);
  int *nelem_part = &_v_nelem_part[0];
  for (int jj=0; jj<nelemfat; jj++) {
    nelem_part[vpart[jj]]++;
    // printf("elem %d in proc %d\n",jj+1,vpart[jj]);
  }

  // Partition nodes. Assign to each node the partition of the first
  // element in tne node2elem list. If there is no element, then
  // assign processor 0 and give a warning.

  // In the future this may be done in a better way. Following links
  // in the other elemensets, even if these elemensets have not been
  // taken into account in the partitioning. 

  int print_nodal_partitioning=0;
  ierr = get_int(mesh->global_options,
		 "print_nodal_partitioning",
		 &print_nodal_partitioning,1); CHKERRQ(ierr);
  if (print_nodal_partitioning) {
#define II_STAT(j,k) VEC2(ii_stat,j,k,size)
    std::vector<int> _v_ii_stat(size*size, 0);
    int *ii_stat = &_v_ii_stat[0];
    PetscPrintf(comm,"\nNodal partitioning (node/processor): \n");
    
    // Node interface between processor statistics
    int c1,c2,P1,P2;
    for (node=0; node<nnod; node++) {
      for (c1 = n2eptr[node]; c1 < n2eptr[node+1]-1; c1++) {
	P1 = vpart[node2elem[c1]];
	for (c2 = c1+1; c2 < n2eptr[node+1]; c2++) {
	  P2 = vpart[node2elem[c2]];
	  if (P1<=P2) II_STAT(P1,P2)++;
	  else II_STAT(P2,P1)++;
	}
      }
    }
    if (size > 1) {
      PetscPrintf(comm,"---\nInter-processor node connections\n");
      for (P1=0; P1<size; P1++)
	for (P2=0; P2<size; P2++) 
	  printf("[%d]-[%d] %d\n",P1,P2,II_STAT(P1,P2));
      PetscPrintf(comm,"\n");
    }
#undef II_STAT
  }

  // node_not_connected_to_fat:= flags whether there is a node not
  // connected to any fat elemset or not. 
  int node_not_connected_to_fat=0;
  int counter=0;
  for (node=0; node<nnod; node++) {
    if (n2eptr[node]==n2eptr[node+1]) {
      node_not_connected_to_fat++;
      npart[node]=1;
    } else {
      int curs_ele = n2eptr[node]; // cursor to element in connected
				   // element list
      int proc;			   // processor to be assigned
#if 1
      // We assign to a node that is connected to elements in
      // different processors the highest processor number. This is
      // done this way for avoiding problems with the IISDMat
      // class. Otherwise an elemset can contribute with local-local
      // elements in other processors, which is not contemplated now. 
      proc = vpart[node2elem[curs_ele]]+1;
      for (curs_ele = n2eptr[node]+1; curs_ele < n2eptr[node+1]; curs_ele++) {
	if (vpart[node2elem[curs_ele]]+1 > proc) 
	  proc = vpart[node2elem[curs_ele]]+1;
      }
#else
      // Take a "random" (connected) processor
      // Number of connected elements
      int conn_ele = n2eptr[node+1] - n2eptr[node];
      int eshift = (counter++) % conn_ele;
      proc = vpart[node2elem[curs_ele + eshift]] +1;
#endif 
      npart[node] = proc;
    }
    if (print_nodal_partitioning)
      PetscPrintf(comm, "%d   %d\n",node+1,npart[node]);
  }
  if (print_nodal_partitioning)
    PetscPrintf(comm,"End nodal partitioning table\n\n");

  if (print_nodal_partitioning && node_not_connected_to_fat)
    PetscPrintf(comm,"warning! there are %d "
		"nodes not linked to any \"fat\" elemset. \n"
		"This induces artificial numbering. [But may be OK]\n",
		node_not_connected_to_fat);

  //o Prints element partitioning. 
  TGETOPTDEF(mesh->global_options,int,debug_element_partitioning,0);

  // Define the eparts of each fat elemset as the corresponding part
  // of vpart. 
  for (std::size_t ielset=0; ielset<nelemsets; ielset++) {
    Elemset* elemset  = *(Elemset **)da_ref(mesh->elemsetlist,ielset);
    if (!elemset->isfat) continue;
    int *epart= elemset->epart;
    for (int iel=0; iel<elemset->nelem; iel++) {
      int ielg=nelemsetptr[ielset]+iel;
      epart[iel] = vpart[ielg]+1;
    }
    if (debug_element_partitioning) {
      PetscPrintf(comm,
		  "Elemset \"%s\", type \"%s\", nelem %d\n."
		  "Partitioning table: \n",elemset->name(),elemset->type,elemset->nelem);
      for (int kk=0; kk<elemset->nelem; kk++) 
	PetscPrintf(comm,"%d -> %d\n",kk,epart[kk]);
    }
  }
  
  delete[] vpart;

  PetscFunctionReturn(0);
}

#undef  __FUNC__
#define __FUNC__ "dofmap_setup"
static int 
dofmap_setup(MPI_Comm comm, Mesh* mesh, Dofmap* dofmap,
	     std::vector<int>& nodepart) {

  PetscErrorCode ierr;
  PetscFunctionBegin;

  assert(mesh   != NULL);
  assert(dofmap != NULL);

  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  int nnod = dofmap->nnod;
  int ndof = dofmap->ndof;
#if 0
  int&   neq       = dofmap->neq;
#else
  int    neq       = dofmap->neq;
#endif
  int*   startproc = dofmap->startproc;
  int*   neqproc   = dofmap->neqproc;
  int*   npart     = &nodepart[0];//dofmap->npart;

  memset(startproc, 0, sizeof(int)*(size+1));
  memset(neqproc,   0, sizeof(int)*(size+1));

  int node, jdof,kdof;
  int nelem, nel, iele, k;
  int *icone;
  row_t row,col; row_t::iterator kndx;

  std::size_t nelemsets = da_length(mesh->elemsetlist);
  int nelemfat=0, is_any_fat=0;
  int *n2eptr, *nelemsetptr;

  // Number degrees of freedom
  int proc;
  std::vector<int> _v_perm(ndof*nnod, 0);
  int *perm = &_v_perm[0];

  // First number all the edof's in processor 0, then on 1, etc...
  // perm contains the permutation. 

  // The criterion to assign a column index (dof or fixed) to a
  // processor is as follows: take the corresponding column, take the
  // first non-null entry, then assign to the processor that the
  // corresponding node belongs. 

  int field;
  map<int,int>::iterator end_fix;
  end_fix = dofmap->fixed_dofs.end();
  
  // First, put perm = to minus the corresponding processor. 
  // The fixed dof's are set to `-(size+1)'  (number of processors + 1)
  for (k=1; k<=nnod*ndof; k++) {
    dofmap->id->get_col(k,col);
    if (col.size()==0) continue;
    if (dofmap->fixed_dofs.find(k) != end_fix) {
      perm[k-1] = -(size+1);
      continue;
    }
    int edof = col.begin()->first;
    dofmap->nodf(edof,node,field);
    perm[k-1] = - npart[node-1];
  }

  // Now, number first all the dofs in processor 0, then on 1, and so
  // on until processor size-1, and finally those fixed (set to perm=size)
  jdof=0;
  for (proc=1; proc<=size+1; proc++) {
    startproc[proc-1]=jdof;
    for (k=1; k<=nnod*ndof; k++) {
      if (perm[k-1] == -proc) perm[k-1] = ++jdof;
    }
    neqproc[proc-1] = jdof-startproc[proc-1];
  }
  neq  = startproc[size];
  dofmap->neq    = neq;
  dofmap->neqf   = neqproc[size];
  dofmap->neqtot = dofmap->neq + dofmap->neqf;
//PetscPrintf(comm,
//	      "Total number of degrees of freedom neq:     %d\n"
//	      "Total number of independent fixations neqf: %d\n",
//	      neq,dofmap->neqf);

//if (size>1) {
//  for (proc=0; proc<size; proc++) 
//    PetscPrintf(comm,
//		  "[%d] number of dof's: %d\n",proc,neqproc[proc]);
//}
  

  idmap *id = new idmap(nnod*ndof,NULL_MAP);
  dofmap->id->remap_cols(perm,*id);
  delete dofmap->id;
  dofmap->id = id;
 
  //o Checks that the  #idmap#  has been correctly generated. 
  TGETOPTDEF(mesh->global_options,int,check_dofmap_id,0);
  if (check_dofmap_id && rank==0) {
    dofmap->id->check();
    dofmap->id->print_statistics();
  }

  //o Prints the dofmap  #idmap#  object. 
  TGETOPTDEF(mesh->global_options,int,print_dofmap_id,0);
  if (print_dofmap_id && rank==0) {
    dofmap->id->print("dofmap->id: \n");
  }

  dofmap->freeze();

  map<int,int>::iterator jj;

  int nfixed = dofmap->fixed.size();
  vector<fixation_entry> fixed_remapped(nfixed);
  //  fixed_remapped.reserve(nfixed);
  for (jj=dofmap->fixed_dofs.begin(); jj!=dofmap->fixed_dofs.end(); jj++) {
    int newjeq = perm[jj->first -1];
    if (newjeq==0) continue;
    int indx= newjeq - neq - 1;
    assert(0<=indx && indx<nfixed); //This should be in this range. Otherwise it falls
				    // off the `fixed_remapped' vector.
    fixed_remapped[indx] =(dofmap->fixed)[jj->second];
  }


  // swap fixed with fixed_remapped (reordered)
  dofmap->fixed.swap(fixed_remapped);
  VOID_IT(fixed_remapped);
  VOID_IT(dofmap->fixed_dofs);

  // Dof's connected to elements in this processor
  int jel,locdof,keq;
  std::vector<int> _dof_here(neq, 0);
  int *dof_here = &_dof_here[0];

  int dof1,dof2; // interval of dof's that live in the processor
  dof1=startproc[rank]+1;
  dof2=dof1+neqproc[rank]-1;
  dofmap->dof1 = dof1;
  dofmap->dof2 = dof2;
  set<int> ghost_dof_set;
    
  // dof_here:= dof_here_list := Auxiliary vectors to define the
  // scatter needed for the ghost values
  for (std::size_t ielset=0; ielset<nelemsets; ielset++) {
    Elemset *elemset  = *(Elemset **)da_ref(mesh->elemsetlist,ielset);
    icone = elemset->icone;
    nelem = elemset->nelem;
    nel = elemset->nel;
    elemset->nelem_here=0;

    // fixme:= Aca hay codigo duplicado. Habria que hacer dos lazos
    // sobre los elemsets. Primero se define el epart para los no-fat
    // y despues se definen los dof_here
    if (elemset->isfat) {
      for (iele=0; iele<nelem; iele++) {
	if(elemset->epart[iele]!=rank+1) continue;
	elemset->nelem_here++;
 	for (jel=0; jel<nel; jel++) {
 	  node = ICONE(iele,jel);
	  for (kdof=1; kdof<=ndof; kdof++) {
	    int m;
	    const int *dofs;
	    const double *coefs;
	    dofmap->get_row(node,kdof,m,&dofs,&coefs);
	    for (int l=0; l<m; l++) {
	      int dof = dofs[l];
	      if (dof <= neq)
		dof_here[dof-1]=1;
	      if (dof <= neq && !(dof1 <= dof && dof <= dof2))
		ghost_dof_set.insert(dof-1);
	    }
	  }
	}
      }
    } else {
      const int *conn; int nell;
      for (iele=0; iele<nelem; iele++) {
	// Decide in which processor will be computed this element
	assert(nel>0);		// Elements should have at least one
				// node
	nell = elemset->real_nodes(iele,conn);
	node = conn[0];
	proc = npart[node-1];
// 	if (proc<1 || proc>size) {
// 	  PetscPrintf(comm,
// 		      "node %d attached to processor %d out of range",node,proc);
// 	  CHKERRA(1);
// 	}
	elemset->epart[iele]=proc;

	if(elemset->epart[iele]!=rank+1) continue;
	elemset->nelem_here++;
	// If decided that belongs to this processor mark all
	// the connected dofs in dof_here
	for (jel=0; jel<nel; jel++) {
	  node = ICONE(iele,jel);
	  for (kdof=1; kdof<=ndof; kdof++) {
	    dofmap->get_row_free(node,kdof,row);
	    for (kndx=row.begin(); kndx!=row.end(); kndx++) {
	      int dof = kndx->first;
	      dof_here[dof-1]=1;
	      if (dof <= neq && !(dof1 <= dof && dof <= dof2))
		ghost_dof_set.insert(dof-1);
	    }
	  }
	}
      }
    }
    { // Make a local scope.
      // Compute the epart2 and epart_p vectors. 
      // Later epart2 could be avoided and only using epart.
      // Right now we leave both, but it is coded so that
      // later we could simply assign epart2=epart and it should work. 
      vector<int> &epart_p = elemset->epart_p;
      epart_p.resize(size+1);
      for (int p=0; p<(size+1); p++) epart_p[p]=0;
      int *epart2 = elemset->epart2;
      int *epart = elemset->epart;
      // First we assign epart2[k] = proc*nelem + <indx-of-k-in-proc>
      // so that we can have the number of processor as epart2[k]/nelem
      for (int k=0; k<nelem; k++) {
	int proc = epart[k]-1;
	assert(proc<size);
	epart2[k] = proc*nelem+epart_p[proc]++;
      }
      int nelem2 = 0;
      for (int p=0; p<size; p++) {
	int tmp = epart_p[p];
	epart_p[p] = nelem2;
	nelem2 += tmp;
      }
      assert(nelem2==nelem);
      epart_p[size] = nelem;
      // Now convert to the true numbering
      for (int k=0; k<nelem; k++) {
	int proc = epart2[k]/nelem;
	epart2[k] += -proc*nelem + epart_p[proc];
      }
      elemset->e1 = epart_p[rank];
      elemset->e2 = epart_p[rank+1];
    }
    // ghost_elems:= These are elements that have related dof's on the
    // processor, but they don't live on the processor.
    // Loops for defining profiles of matrices must loop over them
    // also. 

    if (elemset->ghost_elems == NULL)
      elemset->ghost_elems = da_create(sizeof(int));
    else
      da_resize(elemset->ghost_elems, 0);

    for (iele=0; iele<nelem; iele++) {
      if(elemset->epart[iele]==rank+1) continue;
      for (jel=0; jel<nel; jel++) {
	node = ICONE(iele,jel);
	for (kdof=1; kdof<=ndof; kdof++) {
	  dofmap->get_row(node,kdof,row);
	  for (kndx=row.begin(); kndx!=row.end(); kndx++) {
	    int dof = kndx->first;
	    if (dof1 <= dof && dof <= dof2) {
	      da_append(elemset->ghost_elems,&iele);
	      goto CONTINUE;
	    }
	  }
	}
      }
    CONTINUE:;
    }

    //o Defines a ``locker'' for each element
    TGETOPTDEF(elemset->thash,int,local_store,0);
    if (local_store) {
      if (elemset->local_store) delete[] elemset->local_store;
      elemset->local_store = new void*[elemset->nelem_here];
      for (int j=0; j<elemset->nelem_here; j++) {
	elemset->local_store[j]=NULL;
      }
    }

    da_sort (elemset->ghost_elems,int_cmp,NULL);

  }
  
  // Convert ghost_dof_set en dofmap
  dofmap->ghost_dofs->clear();
  dofmap->ghost_dofs->reserve(ghost_dof_set.size());
  dofmap->ghost_dofs->insert(dofmap->ghost_dofs->end(),
			     ghost_dof_set.begin(),
			     ghost_dof_set.end());
  ghost_dof_set.clear();

  // DEFINE SCATTER
  assert(dofmap->ghost_scatter);
  assert(dofmap->scatter_print);
  if (*dofmap->ghost_scatter != PETSC_NULL) {
    ierr = VecScatterDestroy(*dofmap->ghost_scatter); CHKERRQ(ierr);
    *dofmap->ghost_scatter = PETSC_NULL;
  }
  if (*dofmap->scatter_print != PETSC_NULL) {
    ierr = VecScatterDestroy(*dofmap->scatter_print); CHKERRQ(ierr);
    *dofmap->scatter_print = PETSC_NULL; 
  }
  
  Vec x,ghost_vec,xseq;
  IS is_ghost_glob,is_ghost_loc,is_print;
  PetscInt nghost,*ghosts;
  
  nghost = dofmap->ghost_dofs->size();
  ghosts = &*dofmap->ghost_dofs->begin();
  
  ierr = VecCreateMPI(comm,dofmap->neqproc[rank],neq,&x); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF,nghost,&ghost_vec); CHKERRQ(ierr);

  ierr = ISCreateGeneral(comm,nghost,ghosts,&is_ghost_glob); CHKERRQ(ierr);
  ierr = ISCreateStride(comm,nghost,0,1,&is_ghost_loc);

  ierr = VecScatterCreate(x,is_ghost_glob,ghost_vec,is_ghost_loc,
			  dofmap->ghost_scatter); CHKERRQ(ierr); 

  ierr = ISDestroy(is_ghost_glob); CHKERRQ(ierr); 
  ierr = ISDestroy(is_ghost_loc); CHKERRQ(ierr); 
  ierr = VecDestroy(ghost_vec); CHKERRQ(ierr); 

  int neql = (rank==0) ? dofmap->neq : 0;
  ierr = VecCreateSeq(PETSC_COMM_SELF,neql,&xseq);  CHKERRQ(ierr);
  ierr = ISCreateStride(comm,neql,0,1,&is_print);
  ierr = VecScatterCreate(x,is_print,xseq,is_print,
			  dofmap->scatter_print); CHKERRQ(ierr); 
  
  ierr = ISDestroy(is_print); CHKERRQ(ierr); 

  ierr = VecDestroy(x); CHKERRQ(ierr); 
  ierr = VecDestroy(xseq); CHKERRQ(ierr); 

  PetscFunctionReturn(0);
}


#include "gvars.h"

PF4PY_NAMESPACE_BEGIN

void
setup_mesh(MPI_Comm comm, ::Mesh* mesh, 
	   std::vector<int>&   nodepart,
	   std::vector<float>& weights)
{ 
  assert(mesh != NULL);
  GlobalVars gvars(comm, mesh->global_options, mesh, NULL);
  int ierr = ::mesh_setup(comm, mesh, nodepart, weights);
  if (ierr) throw Error(ierr, "PETSc-FEM error in mesh setup");
}

void
setup_dofmap(MPI_Comm comm, ::Mesh* mesh, ::Dofmap* dofmap,
	     std::vector<int>& nodepart)
{ 
  assert(mesh   != NULL);
  assert(dofmap != NULL);
  GlobalVars gvars(comm, mesh->global_options, mesh, dofmap);
  int ierr = ::dofmap_setup(comm, mesh, dofmap, nodepart); 
  if (ierr) throw Error(ierr, "PETSc-FEM error in dofmap setup");
}

PF4PY_NAMESPACE_END
