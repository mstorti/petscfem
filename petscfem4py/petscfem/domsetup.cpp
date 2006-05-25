//$Id: domsetup.cpp,v 1.1.2.2 2006/05/25 00:34:06 dalcinl Exp $

// STL components
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <cassert>
using namespace std;


#include "petscfem4py.h"
typedef PYPF_NAMESPACE::Error Error;

// PETSc-FEM components
#include <src/fem.h>
#include <src/elemset.h>
#include <src/dofmap.h>

extern Mesh* GLOBAL_MESH;

void metis_part_comm(MPI_Comm comm, int nelemfat,Mesh *mesh,
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
int mesh_setup(MPI_Comm comm, Mesh* mesh, Dofmap* dofmap) {

  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  int nnod = mesh->nodedata->nnod;
  int ndof = dofmap->ndof;

  float* tpwgts    = dofmap->tpwgts;
  int*   npart     = dofmap->npart;

  int ierr;
  int node, jdof,kdof;
  int nelem, nel, iele;
  int *icone;
  row_t row,col; row_t::iterator kndx;

  int nelemsets = da_length(mesh->elemsetlist);
  int nelemfat=0, is_any_fat=0;
  int *n2eptr, *nelemsetptr;

  GLOBAL_MESH = mesh;
  if (size==1) for (int k=0; k<nnod; npart[k++]=1);


  n2eptr = new int[nnod+1];
  nelemsetptr = new int[nelemsets+1];
  nelemsetptr[0]=0;
  for (node=0; node<=nnod; node++) n2eptr[node]=0;

  // compute total number of elements in fat elemsets, nelemfat
  // also n2eptr[node] gets the total number of elements connected to
  // node and n2eptr[0]=0


  for (int i=0; i<nelemsets; i++) {
    Elemset * elemset  = *(Elemset **)da_ref(mesh->elemsetlist,i);
    elemset->ndof = ndof;
    elemset->initialize();
  }

  // numfat:= number of fat elemsets
  for (int ielset=0; ielset<nelemsets; ielset++) {
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
  if (size>1 && !is_any_fat)
    throw Error( "At least one elemset has to be set as fat ('is_fat' option) "
		 "if number of processors is greater than one");

  // cumulated sum 
  for (node=0; node<nnod; node++) 
    n2eptr[node+1] += n2eptr[node];
  // PetscPrintf(comm,"n2eptr[nnod]: %d\n",n2eptr[nnod]);

  // n2esize:= total number of entries in the node2elem table
  int n2esize = n2eptr[nnod];

  // initialize node2elem
  int *node2elem = new int[n2esize];
  for (int j=0; j<n2esize; j++) node2elem[j]=-1;

  // define node2elem pointer
  for (int ielset=0; ielset<nelemsets; ielset++) {
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
  int partflag;
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
    throw Error("partitioning method not known");
  }

  ////PetscPrintf(comm,"Starts partitioning.\n"); 

  int *vpart = new int[nelemfat];
  if (partflag==0 || partflag==2) {

    metis_part_comm(comm, nelemfat,mesh,nelemsets,vpart,
		    nelemsetptr,n2eptr,node2elem,size,rank,
		    partflag,tpwgts,max_partgraph_vertices,
		    iisd_subpart,print_partitioning_statistics);

  } else if (partflag==1) {

    int *mnel = new int[size+1];
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
    delete[] mnel;

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
    ierr = MPI_Bcast (vpart,nelemfat,MPI_INT,0,comm);
  } else assert(0); // something went wrong

//   PetscPrintf(comm,"Ends partitioning.\n");

  // nelem_part:= nelem_part[proc] is the number of elements in
  // processor proc. 
  int *nelem_part = new int[size];
  for (int proc=0; proc<size; proc++) nelem_part[proc]=0;
  for (int jj=0; jj<nelemfat; jj++) {
    nelem_part[vpart[jj]]++;
    // printf("elem %d in proc %d\n",jj+1,vpart[jj]);
  }

  //   if (size>1 && rank==0) {
  //     for (int proc=0; proc<size; proc++) 
  //       printf("%d elements in processor %d\n",nelem_part[proc],proc);
  //   }

  // Partition nodes. Assign to each node the partition of the first
  // element in tne node2elem list. If there is no element, then
  // assign processor 0 and give a warning.

  // In the future this may be done in a better way. Following links
  // in the other elemensets, even if these elemensets have not been
  // taken into account in the partitioning. 

  // node_not_connected_to_fat:= flags whether there is a node not
  // connected to any fat elemset or not. 
  int node_not_connected_to_fat=0;
  
  int print_nodal_partitioning=0;
  ierr = get_int(mesh->global_options,
		 "print_nodal_partitioning",
		 &print_nodal_partitioning,1); CHKERRQ(ierr);
  int counter=0;
#define II_STAT(j,k) VEC2(ii_stat,j,k,size)
  int *ii_stat = new int[size*size];

  if (print_nodal_partitioning)
    PetscPrintf(comm,"\nNodal partitioning (node/processor): \n");

  // Node interface between processor statistics
  int c1,c2,P1,P2;
  for (P1=0; P1<size; P1++) 
    for (P2=0; P2<size; P2++) 
      II_STAT(P1,P2)=0;
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
//   if (rank == 0 && size > 1) {
//     PetscPrintf(comm,"---\nInter-processor node connections\n");
//     for (P1=0; P1<size; P1++) 
//       for (P2=0; P2<size; P2++) 
// 	printf("[%d]-[%d] %d\n",P1,P2,II_STAT(P1,P2));
//     PetscPrintf(comm,"\n");
//   }

  delete[] ii_stat;

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
      PetscPrintf(comm,
		  "%d   %d\n",node+1,npart[node]);
  }
  if (print_nodal_partitioning)
    PetscPrintf(comm,"End nodal partitioning table\n\n");

  if (print_nodal_partitioning && node_not_connected_to_fat)
    PetscPrintf(comm,"warning! there are %d "
		"nodes not linked to any \"fat\" elemset. \n"
		"This induces artificial numbering. [But may be OK]\n",
		node_not_connected_to_fat);

  //o Prints element partitioning. 
  GETOPTDEF(int,debug_element_partitioning,0);

  // Define the eparts of each fat elemset as the corresponding part
  // of vpart. 
  for (int ielset=0; ielset<nelemsets; ielset++) {
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
  
  if (debug_element_partitioning==2) {
    PetscFinalize();
    exit(0);
  }

  delete[] n2eptr;
  delete[] node2elem;
  delete[] nelemsetptr;
  delete[] vpart;
  delete[] nelem_part;
}


#undef  __FUNC__
#define __FUNC__ "dofmap_setup"
int dofmap_setup(MPI_Comm comm, Mesh* mesh, Dofmap* dofmap) {

  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  int nnod = dofmap->nnod;
  int ndof = dofmap->ndof;

  int&   neq       = dofmap->neq;
  int*   startproc = dofmap->startproc;
  int*   neqproc   = dofmap->neqproc;
  float* tpwgts    = dofmap->tpwgts;
  int*   npart     = dofmap->npart;

  int ierr;
  int node, jdof,kdof;
  int nelem, nel, iele, k;
  int *icone;
  row_t row,col; row_t::iterator kndx;

  int nelemsets = da_length(mesh->elemsetlist);
  int nelemfat=0, is_any_fat=0;
  int *n2eptr, *nelemsetptr;

  GLOBAL_MESH = mesh;

  // Number degrees of freedom
  int proc;
  int *perm = new int [ndof*nnod];

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
  dofmap->neq = neq;
  dofmap->neqf = neqproc[size];
  dofmap->neqtot = dofmap->neq + dofmap->neqf;
//   PetscPrintf(comm,
// 	      "Total number of degrees of freedom neq:     %d\n"
// 	      "Total number of independent fixations neqf: %d\n",
// 	      neq,dofmap->neqf);

//   if (size>1) {
//     for (proc=0; proc<size; proc++) 
//       PetscPrintf(comm,
// 		  "[%d] number of dof's: %d\n",proc,neqproc[proc]);
//   }
  

  idmap *idnew = new idmap(nnod*ndof,NULL_MAP);
  dofmap->id->remap_cols(perm,*idnew);
  PYPF_DELETE_SCLR(dofmap->id);
  dofmap->id = idnew;
 
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

  delete[] perm;

  // Dof's connected to elements in this processor
  int jel,locdof,keq;
  int *dof_here = new int [neq];

  for (keq=0; keq<neq; keq++) dof_here[keq]=0;

  int dof1,dof2; // interval of dof's that live in the processor
  dof1=startproc[rank]+1;
  dof2=dof1+neqproc[rank]-1;
  dofmap->dof1 = dof1;
  dofmap->dof2 = dof2;
  set<int> ghost_dof_set;
    
  // dof_here:= dof_here_list := Auxiliary vectors to define the
  // scatter needed for the ghost values
  for (int ielset=0; ielset<nelemsets; ielset++) {
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
      int *epart2 = elemset->epart2;
      int *epart = elemset->epart;
      for (int p=0; p<size+1; p++) epart_p[p]=0;
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

    //elemset->ghost_elems = da_create(sizeof(int));

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
      elemset->local_store = new void*[elemset->nelem_here];
      for (int j=0; j<elemset->nelem_here; j++) {
	elemset->local_store[j]=NULL;
      }
    }

    da_sort (elemset->ghost_elems,int_cmp,NULL);

  }
  
  // Convert ghost_dof_set en dofmap
#if 0
  set<int>::iterator it;
  for (it=ghost_dof_set.begin(); it!=ghost_dof_set.end(); it++) 
    dofmap->ghost_dofs->push_back(*it);
#else
  dofmap->ghost_dofs->clear();
  dofmap->ghost_dofs->reserve(ghost_dof_set.size());
  dofmap->ghost_dofs->insert(dofmap->ghost_dofs->end(),
			     ghost_dof_set.begin(),
			     ghost_dof_set.end());
#endif
  ghost_dof_set.clear();
  delete[] dof_here; 


  // DEFINE SCATTER
  
  Vec x,ghost_vec,xseq;
  IS is_ghost_glob,is_ghost_loc,is_print;
  PetscInt nghost,*ghosts;
  
  PetscFunctionBegin;
  
  nghost = dofmap->ghost_dofs->size();
  ghosts = &*dofmap->ghost_dofs->begin();

  ierr = VecCreateMPI(comm,dofmap->neqproc[rank],neq,&x); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF,nghost,&ghost_vec); CHKERRQ(ierr);

  ierr = ISCreateGeneral(comm,nghost, ghosts, &is_ghost_glob);CHKERRQ(ierr);
  ierr = ISCreateStride(comm,nghost,0,1,&is_ghost_loc);

  ierr = VecScatterCreate(x,is_ghost_glob,ghost_vec,is_ghost_loc,
			  dofmap->ghost_scatter); CHKERRQ(ierr); 

  ierr = ISDestroy(is_ghost_glob); CHKERRQ(ierr); 
  ierr = ISDestroy(is_ghost_loc); CHKERRQ(ierr); 

  int neql = (rank==0 ? dofmap->neq : 0);
  ierr = VecCreateSeq(PETSC_COMM_SELF,neql,&xseq);  CHKERRQ(ierr);
  ierr = ISCreateStride(comm,neql,0,1,&is_print);
  ierr = VecScatterCreate(x,is_print,xseq,is_print,
			  dofmap->scatter_print); CHKERRQ(ierr); 
  
  ierr = ISDestroy(is_print); CHKERRQ(ierr); 

  ierr = VecDestroy(x);CHKERRQ(ierr); 
  ierr = VecDestroy(ghost_vec);CHKERRQ(ierr); 
  ierr = VecDestroy(xseq);CHKERRQ(ierr); 

  PetscFunctionReturn(0);
}


#undef  __FUNC__
#define __FUNC__ "domain_setup"
int domain_setup(MPI_Comm comm, Mesh* mesh, Dofmap* dofmap) {
  return mesh_setup(comm, mesh, dofmap) && dofmap_setup(comm, mesh, dofmap);
}



// #undef  __FUNC__
// #define __FUNC__ "domain_setup"
// int _domain_setup(MPI_Comm comm, Mesh* mesh, Dofmap* dofmap) {

//   int size, rank;
//   MPI_Comm_size(comm, &size);
//   MPI_Comm_rank(comm, &rank);

//   int nnod = dofmap->nnod;
//   int ndof = dofmap->ndof;

//   int&   neq       = dofmap->neq;
//   int*   startproc = dofmap->startproc;
//   int*   neqproc   = dofmap->neqproc;
//   float* tpwgts    = dofmap->tpwgts;
//   int*   npart     = dofmap->npart;

//   int ierr;
//   int node, jdof,kdof;
//   int nelem, nel, iele, k;
//   int *icone;
//   row_t row,col; row_t::iterator kndx;

//   int nelemsets = da_length(mesh->elemsetlist);
//   int nelemfat=0, is_any_fat=0;
//   int *n2eptr, *nelemsetptr;

//   GLOBAL_MESH = mesh;
//   if (size==1) for (k=0; k<nnod; npart[k++]=1);


//   n2eptr = new int[nnod+1];
//   nelemsetptr = new int[nelemsets+1];
//   nelemsetptr[0]=0;
//   for (node=0; node<=nnod; node++) n2eptr[node]=0;

//   // compute total number of elements in fat elemsets, nelemfat
//   // also n2eptr[node] gets the total number of elements connected to
//   // node and n2eptr[0]=0


//   for (int i=0; i<nelemsets; i++) {
//     Elemset * elemset  = *(Elemset **)da_ref(mesh->elemsetlist,i);
//     elemset->ndof = ndof;
//     elemset->initialize();
//   }

//   // numfat:= number of fat elemsets
//   for (int ielset=0; ielset<nelemsets; ielset++) {
//     Elemset * elemset  = *(Elemset **)da_ref(mesh->elemsetlist,ielset);
//     int is_fat=1;
//     ierr = get_int(elemset->thash,"is_fat",&is_fat,1); CHKERRQ(ierr);
//     elemset->isfat=is_fat;
//     nelemsetptr[ielset+1]=nelemsetptr[ielset];
//     if (!is_fat) continue;
//     nelem = elemset->nelem;
//     icone = elemset->icone;
//     nel = elemset->nel;
//     is_any_fat=1;
//     nelemfat += nelem;
//     nelemsetptr[ielset+1]=nelemsetptr[ielset]+nelem;
//     const int *conn; int nell;
//     for (int iel=0; iel<nelem; iel++) {
//       nell = elemset->real_nodes(iel,conn);
//       for (int iloc=0; iloc<nell; iloc++) {
// 	// printf("adding to node %d\n",CONN(iel,iloc));
// 	n2eptr[conn[iloc]]++;
//       }
//     }
//   }

//   // raise error if not defined at least one fat elemset
//   if (size>1 && !is_any_fat)
//     throw Error( "At least one elemset has to be set as fat ('is_fat' option) "
// 		 "if number of processors is greater than one");

//   // cumulated sum 
//   for (node=0; node<nnod; node++) 
//     n2eptr[node+1] += n2eptr[node];
//   // PetscPrintf(comm,"n2eptr[nnod]: %d\n",n2eptr[nnod]);

//   // n2esize:= total number of entries in the node2elem table
//   int n2esize = n2eptr[nnod];

//   // initialize node2elem
//   int *node2elem = new int[n2esize];
//   for (int j=0; j<n2esize; j++) node2elem[j]=-1;

//   // define node2elem pointer
//   for (int ielset=0; ielset<nelemsets; ielset++) {
//     Elemset* elemset  = *(Elemset **)da_ref(mesh->elemsetlist,ielset);
//     if (!elemset->isfat) continue;
//     nelem = elemset->nelem;
//     // icone = elemset->icone;
//     nel = elemset->nel;
//     for (int iel=0; iel<nelem; iel++) {
//       // iel:= global element number cumulated through the elemsetlist,
//       // (only on fat elemsets)
//       int ielg = nelemsetptr[ielset]+iel;
//       const int *conn; int nell;
//       nell = elemset->real_nodes(iel,conn);
//       for (int iloc=0; iloc<nell; iloc++) {
// 	// search for the next free position on the array
// 	node=conn[iloc];
// 	int jpos;
// 	for (jpos=n2eptr[node-1]; jpos<n2eptr[node]; jpos++) 
// 	  if (node2elem[jpos]==-1) break;
// 	if (jpos==n2eptr[node]) {
// 	  printf("node %d, elems %d, range in node2elem %d - %d\n",
// 		 node,n2eptr[node]-n2eptr[node-1],n2eptr[node-1],
// 		 n2eptr[node]-1);
// 	  for (jpos=n2eptr[node-1]; jpos < n2eptr[node]; jpos++)
// 	    printf("%d ",node2elem[jpos]);
// 	  printf("\n");
// 	  assert(0);
// 	}
// 	node2elem[jpos]=ielg;
//       }
//     }
//   }

// // #if 0
// //   for (node=0; node<nnod; node++) {
// //     printf("node %d: ",node+1);
// //     for (int jpos=n2eptr[node]; jpos<n2eptr[node+1]; jpos++) 
// //       printf(" %d",node2elem[jpos]);
// //     printf("\n");
// //   }
// // #endif
//   // partflag = 0 -> METIS partition;
//   // partflag = 1 -> Hitchhiking partition;
//   // partflag = 2 -> Neighbor partition. 
//   // partflag = 3 -> random partition. 
//   // partflag = 4 -> natural partition. 
//   int partflag;
//   //o Set partitioning method. May be set to  #metis# ,
//   //  #hitchhiking# ,  #nearest_neighbor#  or  #random# .
//   TGETOPTDEF_S(mesh->global_options,string,partitioning_method,metis);
//   //o Print graph statistics
//   TGETOPTDEF(mesh->global_options,int,print_partitioning_statistics,0);
// #define INF INT_MAX
//   //o Maximum number of vertices admissible while computing the
//   // partitioning graph.
//   TGETOPTDEF(mesh->global_options,int,max_partgraph_vertices,INF);
// #undef INF
//   // o Number of subpartitions inside each processor. 
//   TGETOPTDEF(mesh->global_options,int,iisd_subpart,1); //nd
//   iisd_subpart = 1; // In order to deactivate subpartitioning at the high level

//   if (partitioning_method == string("metis")) {
//     partflag = 0;
//   } else if (partitioning_method == string("hitchhiking")) {
//     partflag = 1;
//   } else if (partitioning_method == string("nearest_neighbor")) {
//     partflag = 2;
//   } else if (partitioning_method == string("random")) {
//     partflag = 3;
//   } else if (partitioning_method == string("natural")) {
//     partflag = 4;
//   } else {
//     throw Error("partitioning method not known");
//   }

//   ////PetscPrintf(comm,"Starts partitioning.\n"); 

//   int *vpart = new int[nelemfat];
//   if (partflag==0 || partflag==2) {

//     metis_part_comm(comm, nelemfat,mesh,nelemsets,vpart,
// 		    nelemsetptr,n2eptr,node2elem,size,rank,
// 		    partflag,tpwgts,max_partgraph_vertices,
// 		    iisd_subpart,print_partitioning_statistics);

//   } else if (partflag==1) {

//     int *mnel = new int[size+1];
//     if (size>1) {
//       ///if (rank==0) printf("Hitchicking partition, partflag = %d\n",partflag);
//       for (int nproc=0; nproc < size; nproc++) {
// 	mnel[0] = 0;
// 	mnel[nproc+1] = mnel[nproc]+int(nelemfat*tpwgts[nproc]);
// 	if (mnel[size] < nelemfat) mnel[size] = nelemfat-1;
// 	for (int jelem=mnel[nproc]; jelem < mnel[nproc+1]+1; jelem++) { 
// 	  vpart[jelem]=nproc;
// 	}
//       }
//     } else {
//       for (int jjy=0; jjy<nelemfat; jjy++) 
// 	vpart[jjy]=0;
//     }	
//     delete[] mnel;

//   } else if (partflag==3 || partflag==4) {
//     // random partitioning 
//     if (rank==0) {
//       for (int j=0; j < nelemfat; j++) {
// 	vpart[j] = int(drand()*double(size));
// 	if (partflag==4)
// 	  vpart[j] = int(double(j)/double(nelemfat)*double(size));
// 	// Just in case random functions are too close to the limits
// 	// (In theory, rand() should not touch the limits 0, RAND_MAX)
// 	if (vpart[j]>=size) vpart[j]=size-1;
// 	if (vpart[j]<0) vpart[j]=0;
// 	// printf("element %d in [%d]\n",j+1,vpart[j]);
//       }
//     }
//     ierr = MPI_Bcast (vpart,nelemfat,MPI_INT,0,comm);
//   } else assert(0); // something went wrong

// //   PetscPrintf(comm,"Ends partitioning.\n");

//   // nelem_part:= nelem_part[proc] is the number of elements in
//   // processor proc. 
//   int *nelem_part = new int[size];
//   for (int proc=0; proc<size; proc++) nelem_part[proc]=0;
//   for (int jj=0; jj<nelemfat; jj++) {
//     nelem_part[vpart[jj]]++;
//     // printf("elem %d in proc %d\n",jj+1,vpart[jj]);
//   }

//   //   if (size>1 && rank==0) {
//   //     for (int proc=0; proc<size; proc++) 
//   //       printf("%d elements in processor %d\n",nelem_part[proc],proc);
//   //   }

//   // Partition nodes. Assign to each node the partition of the first
//   // element in tne node2elem list. If there is no element, then
//   // assign processor 0 and give a warning.

//   // In the future this may be done in a better way. Following links
//   // in the other elemensets, even if these elemensets have not been
//   // taken into account in the partitioning. 

//   // node_not_connected_to_fat:= flags whether there is a node not
//   // connected to any fat elemset or not. 
//   int node_not_connected_to_fat=0;
  
//   int print_nodal_partitioning=0;
//   ierr = get_int(mesh->global_options,
// 		 "print_nodal_partitioning",
// 		 &print_nodal_partitioning,1); CHKERRQ(ierr);
//   int counter=0;
// #define II_STAT(j,k) VEC2(ii_stat,j,k,size)
//   int *ii_stat = new int[size*size];

//   if (print_nodal_partitioning)
//     PetscPrintf(comm,"\nNodal partitioning (node/processor): \n");

//   // Node interface between processor statistics
//   int c1,c2,P1,P2;
//   for (P1=0; P1<size; P1++) 
//     for (P2=0; P2<size; P2++) 
//       II_STAT(P1,P2)=0;
//   for (node=0; node<nnod; node++) {
//     for (c1 = n2eptr[node]; c1 < n2eptr[node+1]-1; c1++) {
//       P1 = vpart[node2elem[c1]];
//       for (c2 = c1+1; c2 < n2eptr[node+1]; c2++) {
// 	P2 = vpart[node2elem[c2]];
// 	if (P1<=P2) II_STAT(P1,P2)++;
// 	else II_STAT(P2,P1)++;
//       }
//     }
//   }
// //   if (rank == 0 && size > 1) {
// //     PetscPrintf(comm,"---\nInter-processor node connections\n");
// //     for (P1=0; P1<size; P1++) 
// //       for (P2=0; P2<size; P2++) 
// // 	printf("[%d]-[%d] %d\n",P1,P2,II_STAT(P1,P2));
// //     PetscPrintf(comm,"\n");
// //   }

//   delete[] ii_stat;

//   for (node=0; node<nnod; node++) {
//     if (n2eptr[node]==n2eptr[node+1]) {
//       node_not_connected_to_fat++;
//       npart[node]=1;
//     } else {
//       int curs_ele = n2eptr[node]; // cursor to element in connected
// 				   // element list
//       int proc;			   // processor to be assigned
// #if 1
//       // We assign to a node that is connected to elements in
//       // different processors the highest processor number. This is
//       // done this way for avoiding problems with the IISDMat
//       // class. Otherwise an elemset can contribute with local-local
//       // elements in other processors, which is not contemplated now. 
//       proc = vpart[node2elem[curs_ele]]+1;
//       for (curs_ele = n2eptr[node]+1; curs_ele < n2eptr[node+1]; curs_ele++) {
// 	if (vpart[node2elem[curs_ele]]+1 > proc) 
// 	  proc = vpart[node2elem[curs_ele]]+1;
//       }
// #else
//       // Take a "random" (connected) processor
//       // Number of connected elements
//       int conn_ele = n2eptr[node+1] - n2eptr[node];
//       int eshift = (counter++) % conn_ele;
//       proc = vpart[node2elem[curs_ele + eshift]] +1;
// #endif 
//       npart[node] = proc;
//     }
//     if (print_nodal_partitioning)
//       PetscPrintf(comm,
// 		  "%d   %d\n",node+1,npart[node]);
//   }
//   if (print_nodal_partitioning)
//     PetscPrintf(comm,"End nodal partitioning table\n\n");

//   if (print_nodal_partitioning && node_not_connected_to_fat)
//     PetscPrintf(comm,"warning! there are %d "
// 		"nodes not linked to any \"fat\" elemset. \n"
// 		"This induces artificial numbering. [But may be OK]\n",
// 		node_not_connected_to_fat);

//   //o Prints element partitioning. 
//   GETOPTDEF(int,debug_element_partitioning,0);

//   // Define the eparts of each fat elemset as the corresponding part
//   // of vpart. 
//   for (int ielset=0; ielset<nelemsets; ielset++) {
//     Elemset* elemset  = *(Elemset **)da_ref(mesh->elemsetlist,ielset);
//     if (!elemset->isfat) continue;
//     int *epart= elemset->epart;
//     for (int iel=0; iel<elemset->nelem; iel++) {
//       int ielg=nelemsetptr[ielset]+iel;
//       epart[iel] = vpart[ielg]+1;
//     }
//     if (debug_element_partitioning) {
//       PetscPrintf(comm,
// 		  "Elemset \"%s\", type \"%s\", nelem %d\n."
// 		  "Partitioning table: \n",elemset->name(),elemset->type,elemset->nelem);
//       for (int kk=0; kk<elemset->nelem; kk++) 
// 	PetscPrintf(comm,"%d -> %d\n",kk,epart[kk]);
//     }
//   }
  
//   if (debug_element_partitioning==2) {
//     PetscFinalize();
//     exit(0);
//   }

//   delete[] n2eptr;
//   delete[] node2elem;
//   delete[] nelemsetptr;
//   delete[] vpart;
//   delete[] nelem_part;



//   // #####################################################################


//   // Number degrees of freedom
//   int proc;
//   int *perm = new int [ndof*nnod];
//   for (k=0; k<ndof*nnod; k++) perm[k]=0;


//   // First number all the edof's in processor 0, then on 1, etc...
//   // perm contains the permutation. 

//   // The criterion to assign a column index (dof or fixed) to a
//   // processor is as follows: take the corresponding column, take the
//   // first non-null entry, then assign to the processor that the
//   // corresponding node belongs. 

//   int field;
//   map<int,int>::iterator end_fix;
//   end_fix = dofmap->fixed_dofs.end();
  
//   // First, put perm = to minus the corresponding processor. 
//   // The fixed dof's are set to `-(size+1)'  (number of processors + 1)
//   for (k=1; k<=nnod*ndof; k++) {
//     dofmap->id->get_col(k,col);
//     if (col.size()==0) continue;
//     if (dofmap->fixed_dofs.find(k) != end_fix) {
//       perm[k-1] = -(size+1);
//       continue;
//     }
//     int edof = col.begin()->first;
//     dofmap->nodf(edof,node,field);
//     perm[k-1] = - npart[node-1];
//   }

//   // Now, number first all the dofs in processor 0, then on 1, and so
//   // on until processor size-1, and finally those fixed (set to perm=size)
//   jdof=0;
//   for (proc=1; proc<=size+1; proc++) {
//     startproc[proc-1]=jdof;
//     for (k=1; k<=nnod*ndof; k++) {
//       if (perm[k-1] == -proc) perm[k-1] = ++jdof;
//     }
//     neqproc[proc-1] = jdof-startproc[proc-1];
//   }
//   neq  = startproc[size];
//   dofmap->neq = neq;
//   dofmap->neqf = neqproc[size];
//   dofmap->neqtot = dofmap->neq + dofmap->neqf;
// //   PetscPrintf(comm,
// // 	      "Total number of degrees of freedom neq:     %d\n"
// // 	      "Total number of independent fixations neqf: %d\n",
// // 	      neq,dofmap->neqf);

// //   if (size>1) {
// //     for (proc=0; proc<size; proc++) 
// //       PetscPrintf(comm,
// // 		  "[%d] number of dof's: %d\n",proc,neqproc[proc]);
// //   }
  

//   idmap *idnew = new idmap(nnod*ndof,NULL_MAP);
//   dofmap->id->remap_cols(perm,*idnew);
//   PYPF_DELETE_SCLR(dofmap->id);
//   dofmap->id = idnew;
 
//   //o Checks that the  #idmap#  has been correctly generated. 
//   TGETOPTDEF(mesh->global_options,int,check_dofmap_id,0);
//   if (check_dofmap_id && rank==0) {
//     dofmap->id->check();
//     dofmap->id->print_statistics();
//   }

//   //o Prints the dofmap  #idmap#  object. 
//   TGETOPTDEF(mesh->global_options,int,print_dofmap_id,0);
//   if (print_dofmap_id && rank==0) {
//     dofmap->id->print("dofmap->id: \n");
//   }

//   dofmap->freeze();

//   map<int,int>::iterator jj;

// // #if 0  // debug:=
// //   {
// //     map<int,int>::iterator jjjj;
// //     for (jjjj=dofmap->fixed_dofs.begin(); jjjj!=end_fix; jjjj++) 
// //       printf("%d -> fix %d\n",jjjj->first,jjjj->second);

// //     for (int kkk=0; kkk<nnod*ndof; kkk++) {
// //       printf("%d -> perm[%d]\n",kkk,perm[kkk]);
// //     }
  
// //   }
// // #endif

//   int nfixed = dofmap->fixed.size();
//   vector<fixation_entry> fixed_remapped(nfixed);
//   //  fixed_remapped.reserve(nfixed);
//   for (jj=dofmap->fixed_dofs.begin(); jj!=dofmap->fixed_dofs.end(); jj++) {
//     int newjeq = perm[jj->first -1];
//     if (newjeq==0) continue;
//     int indx= newjeq - neq - 1;
//     assert(0<=indx && indx<nfixed); //This should be in this range. Otherwise it falls
// 				    // off the `fixed_remapped' vector.
//     fixed_remapped[indx] =(dofmap->fixed)[jj->second];
//   }


//   // swap fixed with fixed_remapped (reordered)
//   dofmap->fixed.swap(fixed_remapped);
//   VOID_IT(fixed_remapped);
//   VOID_IT(dofmap->fixed_dofs);

// // #if 0
// //   //debug:= Verificar como quedaron las fijaciones
// //   for (int j=dofmap->neq+1; j<=dofmap->neqtot; j++) {
// //     dofmap->id->get_col(j,col);
// //     assert(col.size()==1);
// //     edof = col.begin()->first;
// //     assert(col.begin()->second == 1.);
// //     dofmap->nodf(edof,node,field);
// //     printf("node %d, field %d , edof %d, j %d, val %f\n",
// // 	   node, field , edof,j,fixed_remapped[j-neq-1]);
// //   }
// // #endif
//   delete[] perm;

//   // Dof's connected to elements in this processor
//   int jel,locdof,keq;
//   int *dof_here = new int [neq];

//   for (keq=0; keq<neq; keq++) dof_here[keq]=0;

//   int dof1,dof2; // interval of dof's that live in the processor
//   dof1=startproc[rank]+1;
//   dof2=dof1+neqproc[rank]-1;
//   dofmap->dof1 = dof1;
//   dofmap->dof2 = dof2;
//   set<int> ghost_dof_set;
    
//   // dof_here:= dof_here_list := Auxiliary vectors to define the
//   // scatter needed for the ghost values
//   for (int ielset=0; ielset<nelemsets; ielset++) {
//     Elemset *elemset  = *(Elemset **)da_ref(mesh->elemsetlist,ielset);
//     icone = elemset->icone;
//     nelem = elemset->nelem;
//     nel = elemset->nel;
//     elemset->nelem_here=0;

//     // fixme:= Aca hay codigo duplicado. Habria que hacer dos lazos
//     // sobre los elemsets. Primero se define el epart para los no-fat
//     // y despues se definen los dof_here
//     if (elemset->isfat) {
//       for (iele=0; iele<nelem; iele++) {
// 	if(elemset->epart[iele]!=rank+1) continue;
// 	elemset->nelem_here++;
//  	for (jel=0; jel<nel; jel++) {
//  	  node = ICONE(iele,jel);
// 	  for (kdof=1; kdof<=ndof; kdof++) {
// 	    int m;
// 	    const int *dofs;
// 	    const double *coefs;
// 	    dofmap->get_row(node,kdof,m,&dofs,&coefs);
// 	    for (int l=0; l<m; l++) {
// 	      int dof = dofs[l];
// 	      if (dof <= neq)
// 		dof_here[dof-1]=1;
// 	      if (dof <= neq && !(dof1 <= dof && dof <= dof2))
// 		ghost_dof_set.insert(dof-1);
// 	    }
// 	  }
// 	}
//       }
//     } else {
//       const int *conn; int nell;
//       for (iele=0; iele<nelem; iele++) {
// 	// Decide in which processor will be computed this element
// 	assert(nel>0);		// Elements should have at least one
// 				// node
// 	nell = elemset->real_nodes(iele,conn);
// 	node = conn[0];
// 	proc = npart[node-1];
// // 	if (proc<1 || proc>size) {
// // 	  PetscPrintf(comm,
// // 		      "node %d attached to processor %d out of range",node,proc);
// // 	  CHKERRA(1);
// // 	}
// 	elemset->epart[iele]=proc;

// 	if(elemset->epart[iele]!=rank+1) continue;
// 	elemset->nelem_here++;
// 	// If decided that belongs to this processor mark all
// 	// the connected dofs in dof_here
// 	for (jel=0; jel<nel; jel++) {
// 	  node = ICONE(iele,jel);
// 	  for (kdof=1; kdof<=ndof; kdof++) {
// 	    dofmap->get_row_free(node,kdof,row);
// 	    for (kndx=row.begin(); kndx!=row.end(); kndx++) {
// 	      int dof = kndx->first;
// 	      dof_here[dof-1]=1;
// 	      if (dof <= neq && !(dof1 <= dof && dof <= dof2))
// 		ghost_dof_set.insert(dof-1);
// 	    }
// 	  }
// 	}
//       }
//     }
//     { // Make a local scope.
//       // Compute the epart2 and epart_p vectors. 
//       // Later epart2 could be avoided and only using epart.
//       // Right now we leave both, but it is coded so that
//       // later we could simply assign epart2=epart and it should work. 
//       vector<int> &epart_p = elemset->epart_p;
//       epart_p.resize(size+1);
//       int *epart2 = elemset->epart2;
//       int *epart = elemset->epart;
//       for (int p=0; p<size+1; p++) epart_p[p]=0;
//       // First we assign epart2[k] = proc*nelem + <indx-of-k-in-proc>
//       // so that we can have the number of processor as epart2[k]/nelem
//       for (int k=0; k<nelem; k++) {
// 	int proc = epart[k]-1;
// 	assert(proc<size);
// 	epart2[k] = proc*nelem+epart_p[proc]++;
//       }
//       int nelem2 = 0;
//       for (int p=0; p<size; p++) {
// 	int tmp = epart_p[p];
// 	epart_p[p] = nelem2;
// 	nelem2 += tmp;
//       }
//       assert(nelem2==nelem);
//       epart_p[size] = nelem;
//       // Now convert to the true numbering
//       for (int k=0; k<nelem; k++) {
// 	int proc = epart2[k]/nelem;
// 	epart2[k] += -proc*nelem + epart_p[proc];
//       }
//       elemset->e1 = epart_p[rank];
//       elemset->e2 = epart_p[rank+1];
// // #if 0
// //       if (!rank) {
// // 	for (int k=0; k<elemset->nelem; k++) {
// // 	  printf("elem %d, epart %d, epart2 %d\n",k,epart[k],epart2[k]);
// // 	}
// //       }
// //       PetscFinalize();
// //       exit(0);
// // #endif    
//     }
//     // ghost_elems:= These are elements that have related dof's on the
//     // processor, but they don't live on the processor.
//     // Loops for defining profiles of matrices must loop over them
//     // also. 

//     //elemset->ghost_elems = da_create(sizeof(int));

//     for (iele=0; iele<nelem; iele++) {
//       if(elemset->epart[iele]==rank+1) continue;
//       for (jel=0; jel<nel; jel++) {
// 	node = ICONE(iele,jel);
// 	for (kdof=1; kdof<=ndof; kdof++) {
// 	  dofmap->get_row(node,kdof,row);
// 	  for (kndx=row.begin(); kndx!=row.end(); kndx++) {
// 	    int dof = kndx->first;
// 	    if (dof1 <= dof && dof <= dof2) {
// 	      da_append(elemset->ghost_elems,&iele);
// 	      goto CONTINUE;
// 	    }
// 	  }
// 	}
//       }
//     CONTINUE:;
//     }

//     //o Defines a ``locker'' for each element
//     TGETOPTDEF(elemset->thash,int,local_store,0);
//     if (local_store) {
//       elemset->local_store = new void*[elemset->nelem_here];
//       for (int j=0; j<elemset->nelem_here; j++) {
// 	elemset->local_store[j]=NULL;
//       }
//     }

//     da_sort (elemset->ghost_elems,int_cmp,NULL);
//     int nghostel = da_length(elemset->ghost_elems);

// //     if (size>1) {
// //       PetscPrintf(comm,
// // 		  "For elemset type \"%s\", name \"%s\"\n",
// // 		  elemset->type,elemset->name());
// //       PetscSynchronizedPrintf(comm,
// // 			      "On processor [%d], %d local elements,"
// // 			      " %d ghost elements.\n",
// // 			      rank,elemset->nelem_here,nghostel);
// //       PetscSynchronizedFlush(comm);
// //     }

//   }
  
//   // Convert ghost_dof_set en dofmap
//   set<int>::iterator it;
//   for (it=ghost_dof_set.begin(); it!=ghost_dof_set.end(); it++) 
//     dofmap->ghost_dofs->push_back(*it);
//   PetscSynchronizedFlush(comm);
//   VOID_IT(ghost_dof_set);
//   int nghost_dofs = dofmap->ghost_dofs->size();

//   // This dof_here stuff may be done with STL sets
//   int ndofhere = 0;
//   for (keq=0; keq<neq; keq++) {
//     ndofhere += dof_here[keq];
//   }

//   int *dof_here_list;
//   dof_here_list = new int[ndofhere];

//   k=0;
//   for (keq=0; keq<neq; keq++) {
//     if (dof_here[keq]) {
//       dof_here_list[k++] = keq;
//     }
//   }
//   // fixme:= Al final ya practicamente saque todo lo de dof_here y
//   // dof_here_list, asi que esto de arriba y bajo habria que sacarlo. 
//   delete[] dof_here; 
//   delete[] dof_here_list;

//   // Defines certain quantities in dofmap
//   //   dofmap->neq= neq;
//   //   dofmap->startproc = startproc;
//   //   dofmap->neqproc = neqproc;
//   //   dofmap->size = size;
//   //   dofmap->npart = npart;
  
//   //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
//   //                        DEFINE SCATTER
//   //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

// //   PetscPrintf(comm,"Defining scatters...\n");

//   Vec x,ghost_vec,xseq;
//   IS is_ghost_glob,is_ghost_loc,is_print;

//   dofmap->create_MPI_vector(x);
//   dofmap->create_MPI_ghost_vector(ghost_vec);

//   int nghost_tot;
//   ierr = VecGetSize(ghost_vec,&nghost_tot); CHKERRQ(ierr);

//   ierr = ISCreateGeneral(comm,nghost_dofs,
// 			 &*dofmap->ghost_dofs->begin(),
// 			 &is_ghost_glob);  CHKERRQ(ierr); 
//   ierr = ISCreateStride(comm,nghost_dofs,0,1,&is_ghost_loc);

//   ierr = VecScatterCreate(x,is_ghost_glob,ghost_vec,is_ghost_loc,
// 			  dofmap->ghost_scatter); CHKERRQ(ierr); 

//   ierr = ISDestroy(is_ghost_glob); CHKERRQ(ierr); 
//   ierr = ISDestroy(is_ghost_loc); CHKERRQ(ierr); 

//   int neql = (rank==0 ? dofmap->neq : 0);
//   ierr = VecCreateSeq(PETSC_COMM_SELF,neql,&xseq);  CHKERRQ(ierr);
  
//   ierr = ISCreateStride(comm,neql,0,1,&is_print);
//   ierr = VecScatterCreate(x,is_print,xseq,is_print,
// 			  dofmap->scatter_print); CHKERRQ(ierr); 
  
//   ierr = ISDestroy(is_print); CHKERRQ(ierr); 

// // #if 0
// //   for (int jj=dof1; jj<=dof2; jj++) {
// //     VecSetValue(x,jj-1,double(jj),INSERT_VALUES);
// //   }
// //   PetscPrintf(comm,"vector x en read_mesh\n");
// //   ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD); CHKERRA(ierr);
// //   PetscPrintf(comm,"================\n");
// //   ierr = VecScatterBegin(x,ghost_vec,INSERT_VALUES,
// // 			 SCATTER_FORWARD,dofmap->ghost_scatter); CHKERRA(ierr); 
// //   ierr = VecScatterEnd(x,ghost_vec,INSERT_VALUES,
// // 		       SCATTER_FORWARD,dofmap->ghost_scatter); CHKERRA(ierr); 
  
// //   double *array;
// //   ierr = VecGetArray(ghost_vec,&array); CHKERRQ(ierr);

// //   PetscSynchronizedPrintf(comm,"Local ghost values on [%d]\n",rank);
// //   for (int jj=0; jj<nghost_dofs; jj++ ) 
// //     PetscSynchronizedPrintf(comm,"local %d, dof %d  -> %g\n",
// // 			    jj,(*dofmap->ghost_dofs)[jj],array[jj]);
// //   PetscSynchronizedFlush(comm);
// // #endif

//   ierr = VecDestroy(x);
//   ierr = VecDestroy(ghost_vec);
//   ierr = VecDestroy(xseq);

// // #if 0
// //   ierr = VecGhostGetLocalForm(gx,&lx); CHKERRQ(ierr);

// //   for (int jj=dof1; jj<=dof2; jj++) {
// //     VecSetValue(gx,jj-1,double(jj),INSERT_VALUES);
// //   }
// //   ierr = VecAssemblyBegin(gx); CHKERRQ(ierr);
// //   ierr = VecAssemblyEnd(gx); CHKERRQ(ierr);

// //   ierr = VecGhostUpdateBegin(gx,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
// //   ierr = VecGhostUpdateEnd(gx,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

// //   double *array;
// //   ierr = VecGetArray(lx,&array); CHKERRQ(ierr);

// //   PetscSynchronizedPrintf(comm,"On processor [%d], local values\n",rank);
// //   for (int jj=0; jj<neqproc[rank]; jj++ ) {
// //     PetscSynchronizedPrintf(comm,"local %d, global %d -> %g\n",
// // 			    jj,startproc[rank]+jj,array[jj]);
// //   }

// //   PetscSynchronizedPrintf(comm,"Local ghost values\n");
// //   for (int jj=0; jj<nghost_dofs; jj++ ) {
// //     int local = neqproc[rank]+jj;
// //     PetscSynchronizedPrintf(comm,"local %d, global %d -> %g\n",
// // 			    local,(*(dofmap->ghost_dofs))[jj],array[local]);
// //   }

// //   PetscSynchronizedFlush(comm);
// //   ierr = VecRestoreArray(lx,&array);CHKERRQ(ierr);
// //   ierr = VecGhostRestoreLocalForm(gx,&lx);CHKERRQ(ierr); 
// //   PetscFinalize();
// //   exit(0);
// // #endif 

// // #if 0
// //   // para debuggear
// //   ierr = VecScatterBegin(x,xseq,INSERT_VALUES,SCATTER_FORWARD,scatter); CHKERRQ(ierr); 
// //   ierr = VecScatterEnd(x,xseq,INSERT_VALUES,SCATTER_FORWARD,scatter); CHKERRQ(ierr); 
// //   PetscPrintf(comm,"Despues del scatter\n");


// //   // debug
// //   ierr = VecGetArray(xseq,&sstate); CHKERRQ(ierr); 
// //   PetscSynchronizedPrintf(comm,"On processor %d ------------------\n",
// // 			  rank+1);
// //   for (int k=0; k<neq; k++) {
// //     PetscSynchronizedPrintf(comm,"%d -> %f\n",k+1,sstate[k]);
// //   }
// //   PetscSynchronizedFlush(comm);
// //   ierr = VecRestoreArray(xseq,&sstate); CHKERRQ(ierr); 

// //   PetscFinalize();
// //   exit(0);
// // #endif  

// // #if 0
// //   PetscSynchronizedPrintf(comm,
// // 			  "On processor %d \n",rank+1);
// //   for (k=0; k<ndofhere; k++) {
// //     PetscSynchronizedPrintf(comm,
// // 			    "%d -> %d\n",k+1,dof_here_list[k]);
// //   }
// //   PetscSynchronizedFlush(comm);
// // #endif
  
// //   PetscPrintf(comm,"Total number of dof's: %d\n",neq);

// // #if 0
// //   int nloads=0;
// //   fstack->get_line(line);
// //   if (strncmp("loads",line,5) ) {
// //     PetscPrintf(comm,"Couldn't find <fixa> tag line\n");
// //     exit(1);
// //   }
// //   while (1) {
// //     fstack->get_line(line);
// //     if (strstr(line,"__END_LOADS__")) break;
// //     nloads++;
// //     int nr=sscanf(line,"%d %d %lf",&node,&kdof,&dval);
// //     if (nr !=3) {
// //       PetscPrintf(comm,"Error reading LOADS, for loads %d\n",nfixa);
// //     }
// //     //LOAD(node-1,kdof)=dval;
// //     CHKERRQ(1);
// //   }
// //   PetscPrintf(comm,"Total loads: %d\n",nloads);
// // #endif

//   return 0;
// }
