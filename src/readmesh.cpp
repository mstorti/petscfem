//__INSERT_LICENSE__
//$Id: readmesh.cpp,v 1.17 2001/05/05 13:13:26 mstorti Exp $
 
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
#include <algorithm>
#include <cassert>

#include "getprop.h"

using namespace std;

#define IDENT(j,k) VEC2(ident,j,k,ndof)

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
#undef ICONE
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define ELEMIPROPS(j,k) VEC2(elemiprops,j,k,neliprops)
#define NODEDATA(j,k) VEC2(mesh->nodedata->nodedata,j,k,nu)

// fixme:= aca no se porque tuve que pasar neq por referenciar
// porque sino no pasaba correctamente el valor. 
#undef __FUNC__
#define __FUNC__ "read_mesh"
int read_mesh(Mesh *& mesh,char *fcase,Dofmap *& dofmap,
	      int & neq,int size,int myrank) {

  vector<Amplitude *> amplitude_list;

  char *p1,*p2, *token, *type, *bsp=" \t";
  char *line;
  const char *cline;
  Autostr *linecopy = astr_create();
  char *key,*val;
  int ndim,nu,ndof,nnod,ierr,numfat,node,jdof,kdof,edof;
  int pos, nelem, nel, nelprops, neliprops, nread, elemsetnum=0,
	fat_flag,iele,k,nfixa, *ident,rflag;
  double *dptr,dval; 
  TextHashTable *thash;
  Nodedata *nodedata;
  // este despues habria que borrarlo
  map<int,int> fixed_dofs;

  row_t row,col,row0;
  row_t::iterator kndx;

  Elemset *elemset;
  int *icone,etype,*dof_here;
  mesh = new Mesh;
  mesh->nodedata = new Nodedata;
  mesh->elemsetlist = da_create(sizeof(Elemset *));

  // Read data
  FileStack *fstack;
  fstack = new FileStack(fcase);
  if (myrank==0) fstack->set_echo_stream(stdout);

  // Initialize number of eqs.
  neq=0;
  dofmap = new Dofmap;

  while (1) {
    if(fstack->get_line(line)) break;
      
    token = strtok(line,bsp);

    if (token==NULL) {

      continue;

    } else if (!strcmp(token,"global_options")) {
  
      // Reads general data hash table
      read_hash_table(fstack,mesh->global_options);
      if (myrank==0) mesh->global_options->print("\n -- Global_options: ");
      mesh->global_options->register_name("global_options");
      mesh->global_options->set_as_global();

    } else if (!strcmp(token,"nodes")) {

      PetscPrintf(PETSC_COMM_WORLD," -- Reading nodes:\n");
      token = strtok(NULL,bsp); sscanf(token,"%d",&ndim);
      token = strtok(NULL,bsp); sscanf(token,"%d",&nu);
      token = strtok(NULL,bsp); sscanf(token,"%d",&ndof);
      PetscPrintf(PETSC_COMM_WORLD, 
		  "Dimension: %d, Size of nodedata vector: %d\n",ndim,nu);
      mesh->nodedata->nu = nu;
      mesh->nodedata->ndim = ndim;

      dofmap->ndof = ndof;
      node = 0;
      double *row = new double[nu];
      darray *xnod;
      xnod = da_create(nu*sizeof(double));
      while (1) {
	fstack->get_line(line);
	if (strstr("__END_NODES__",line)) break;
	node++;
	for (int kk=0; kk<nu; kk++) {
	  token = strtok((kk==0 ? line : NULL),bsp);
	  if (token==NULL) {
	    PetscPrintf(PETSC_COMM_WORLD,
			"Error reading coordinates in line:\n\"%s\"\n"
			"Not enough values in line!!\n",line);
	    CHKERRQ(1);
	  }
	  int nread = sscanf(token,"%lf",row+kk);
	  if (nread != 1) {
	    PetscPrintf(PETSC_COMM_WORLD,
			"Error reading coordinates in line:\n\"%s\"",line);
	    CHKERRQ(1);
	  }
	}
	int indx = da_append (xnod,row);
	if (indx<0) PFEMERRQ("Insufficient memory reading nodes");
      }

      // wait_from_console("after reading nodes"); 

      nnod=node;
      dofmap->nnod = nnod;
      mesh->nodedata->nnod = nnod;
      PetscPrintf(PETSC_COMM_WORLD,"Read %d nodes\n",nnod);
      
      delete[] row;
      mesh->nodedata->nodedata = new double[nnod*nu];
      for (node=0; node<nnod; node++) {
	row = (double *) da_ref(xnod,node);
	for (int kk=0; kk<nu; kk++) {
	  NODEDATA(node,kk) = row[kk];
	}
      }
      da_destroy(xnod);

      // calling dofmap constructor
      dofmap->id = new idmap(nnod*ndof,NULL_MAP);

    } else if (!strcmp(token,"table")) {

      token = strtok(NULL,bsp);
      string name(token);
      read_hash_table(fstack,thash);
      thash->register_name(name);
      if (myrank==0) {
	printf("\n\n -- Table \"%s\" read:",name.c_str());
	thash->print("");
      }
      if (name=="global_options") thash->set_as_global();
      mesh->global_options = thash;

    } else if (!strcmp(token,"elemset")) {

      // Now read elemset's
      fat_flag=0;
      token = strtok(NULL,bsp);
      type = new char[strlen(token)+1];
      strcpy(type,token);
      PetscPrintf(PETSC_COMM_WORLD,"\n -- Reading elemset type \"%s\"\n",type);
      sscanf(strtok(NULL,bsp),"%d",&nel);
      token = strtok(NULL,bsp);

      // Reads hash table for the elemeset
      read_hash_table(fstack,thash);
      const char *name;
      thash->get_entry("name",name);
      if (name) {
	thash->register_name(name);
      } else {
	PetscPrintf(PETSC_COMM_WORLD,
		    "No name for this elemset!! setting as \"__ANONYMOUS__\"...\n");
	thash->register_name("__ANONYMOUS__");
      }
      if (myrank==0) thash->print("Table of properties:");

      // Read props line
      char *buf;
      vector<string> prop_names;
      vector<int> prop_lens;
      thash->get_entry("props",cline);
      GHashTable *props =
	g_hash_table_new(&g_str_hash,&g_str_equal);
      nelprops=0;
      int posit=0;
      const char *prop_name;
      if (cline!=NULL) {
	parse_props_line(cline,prop_names,prop_lens);
	props_hash_entry *phe, *pheold;
	for (int j=0; j<prop_names.size(); j++) {
	  prop_name = prop_names[j].c_str();
	  phe = new props_hash_entry;
	  // fixme:= despues hay que poner width a la cantidad de
	  // escalares que contiene esta propiedad
	  phe->width=prop_lens[j];
	  phe->position = posit;
	  posit += phe->width;

	  // verificar que la propiedad no este duplicada
	  pheold = (props_hash_entry *)g_hash_table_lookup(props,prop_name);
	  if (pheold) {
	    PetscPrintf(PETSC_COMM_WORLD,
			"duplicated elem properties label \"%s\"\n",prop_name);
	    PetscFinalize();
	    exit(1);
	  }
	  key = new char[strlen(prop_name)+1];
	  strcpy(key,prop_name);
	  g_hash_table_insert (props,(void *)key,phe);
	}
	g_hash_table_freeze(props);
      }
      nelprops = posit;

      // Read iprops line
      thash->get_entry("iprops",cline);
      neliprops=0; 
      GHashTable *iprops = g_hash_table_new(&g_str_hash,&g_str_equal);
      // GHashTable *iprops = g_hash_table_new(&hash_func,&key_compare_func);
      if (cline!=NULL) {
	//astr_copy_s(linecopy,line);
	buf = local_copy(cline);
	props_hash_entry *phe, *pheold;
	int posit=0;
	while (1) {
	  token = strtok((neliprops==0? buf : NULL),bsp);
	  if (token==NULL) break;
	  neliprops++;
	  phe = new props_hash_entry;
	  // fixme:= despues hay que poner width a la cantidad de
	  // escalares que contiene esta propiedad
	  phe->width=1;
	  phe->position = posit;
	  posit += phe->width;

	  // verificar que la propiedad no este duplicada
	  pheold = (props_hash_entry *)g_hash_table_lookup(props,token);
	  if (pheold) {
	    PetscPrintf(PETSC_COMM_WORLD,
			"duplicated elem properties label \"%s\"\n",token);
	    PetscFinalize();
	    exit(1);
	  }
	  key = new char[strlen(token)+1];
	  strcpy(key,token);
	  g_hash_table_insert (props,(void *)key,phe);
	}
	g_hash_table_freeze(props);
	delete[] buf;
      }

      // read connectivity and element properties
      Autobuf *buff;
      Darray *da_icone;
      int *icorow = new int[nel];
      double *proprow = new double[nelprops];
      int *iproprow = new int[neliprops];
      int rowsize = (nel+neliprops)*sizeof(int)+nelprops*sizeof(double);
      da_icone = da_create(rowsize);
      buff = abuf_create();
      iele=0;
      int node;
      while (1) {
	fstack->get_line(line);

	// reading element connectivities
	for (int jel=0; jel<nel; jel++) {
	  token =  strtok(( jel==0 ? line : NULL),bsp);
	  if (jel==0 && !strcmp(token,"__END_ELEMSET__"))
	    goto DONE;
	  sscanf(token ,"%d",&node);
	  icorow[jel]= node;
	  // This should be done AFTER reading the nodes 
	  // Set all nodes that are connected to an element as degrees of freedom
	  for (int kdof=1; kdof<=ndof; kdof++) {
	    edof = dofmap->edof(node,kdof);
	    dofmap->id->set_elem(edof,edof,1.);
	  }
	}

	// reading element properties
	for (int jprop=0; jprop<nelprops; jprop++) {
	  token = strtok(NULL,bsp);
	  if (token==NULL) {
	    PetscPrintf(PETSC_COMM_WORLD,
			"fails to read per-element property %d,\n at line \"%s\"",
		   jprop+1,line);
	    PFEMERRQ("");
	  }
	  sscanf(token,"%lf",proprow+jprop);
	}

	// reading integer element properties
	for (int jprop=0; jprop<neliprops; jprop++) {
	  token = strtok(NULL,bsp);
	  if (token==NULL) {
	    PetscPrintf(PETSC_COMM_WORLD,
			"fails to read integer per-element"
			" property %d,\n at line \"%s\"",
		   jprop+1,line);
	    PFEMERRQ("");
	  }
	  sscanf(token,"%d",iproprow+jprop);
	}
	
	// Copying to buffer
	abuf_zero (buff);
	abuf_cat_buf (buff,(unsigned char *)icorow,nel*sizeof(int));
	abuf_cat_buf (buff,(unsigned char *)proprow,nelprops*sizeof(double));
	abuf_cat_buf (buff,(unsigned char *)iproprow,neliprops*sizeof(int));
	unsigned char *pp = abuf_data(buff);

	int indxi = da_append(da_icone,abuf_data(buff));
	if ( indxi==-1 ) PFEMERRQ("Insufficient memory reading elements");
      }

    DONE:
      elemsetnum++;
      delete[] icorow;

      nelem = da_length(da_icone);

      double *elemprops = new double[nelem*nelprops];
      int *elemiprops = new int[nelem*neliprops];
      icone = new int[nel*nelem];
      unsigned char *buffp;
      
      for (iele=0; iele<nelem; iele++) {
	icorow=(int *) da_ref(da_icone,iele);
	for (int kk=0; kk<nel; kk++) {
	  ICONE(iele,kk) = icorow[kk];
	}

	buffp = (unsigned char *)da_ref(da_icone,iele)+nel*sizeof(int);
	memcpy ((void *)proprow,buffp,nelprops*sizeof(double));
	for (int kk=0; kk<nelprops; kk++) {
	  ELEMPROPS(iele,kk) = proprow[kk];
	}

	buffp += nelprops*sizeof(double);
	memcpy ((void *)iproprow,buffp,neliprops*sizeof(int));
	for (int kk=0; kk<neliprops; kk++) {
	  ELEMIPROPS(iele,kk) = iproprow[kk];
	}

      }

      delete[] proprow;
      delete[] iproprow;
      da_destroy(da_icone);
      abuf_destroy(buff);

      int nelprops_add, neliprops_add;
      double *elemprops_add;
      int *elemiprops_add;
      //o Additional properties (used by the element routine)
      TGETOPTDEF(thash,int,additional_props,0);
      nelprops_add = additional_props;
      if (nelprops_add>0) {
	elemprops_add = new
	  double[nelem*nelprops_add];
	for (int k=0; k<nelem*nelprops_add; k++)
	  elemprops_add[k]=0.;
      }

      //o int additional properties (used by the element routine)
      TGETOPTDEF(thash,int,additional_iprops,0);
      neliprops_add = additional_iprops;
      if (neliprops_add>0) {
	elemiprops_add = new 
	  int[nelem*neliprops_add];
	for (int k=0; k<nelem*neliprops_add; k++)
	  elemiprops_add[k]=0;
      }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
      // Bless with the appropriate type
      // adding a line CHECK_ELEMSET_TYPE(type) for each
      // type of elemset to be added
      bless_elemset(type,elemset);
      PetscPrintf(PETSC_COMM_WORLD,
		  "elemset number %d, pointer %p, number of elements %d\n",
		  elemsetnum,elemset,nelem);
      
      elemset->type = type;
      elemset->nelem = nelem; 
      elemset->nel   = nel  ; 
      elemset->elem_conne = new int[nel];
      elemset->ndof  = ndof ; 
      elemset->nelprops = nelprops; 
      elemset->neliprops = neliprops; 
      elemset->nelprops_add = nelprops_add; 
      elemset->neliprops_add = neliprops_add; 
      elemset->elemprops_add = elemprops_add; 
      elemset->elemiprops_add = elemiprops_add; 
      elemset->thash = thash; 
      elemset->icone = icone; 
      elemset->elemprops  = elemprops; 
      elemset->elemiprops  = elemiprops; 
      elemset->elem_prop_names  = props; 
      elemset->epart = NULL;
      elemset->isfat = 0;

      // Append to the list
      da_append(mesh->elemsetlist,&elemset);
      PetscPrintf(PETSC_COMM_WORLD,"Ends reading  elemset\n");

    } else if (!strcmp(token,"end_elemsets")) {

      // nothing is done here
      PetscPrintf(PETSC_COMM_WORLD,"End elemsets section\n");
      // wait_from_console("reading end_elemset"); 

    } else if (!strcmp(token,"fixa")) {

      PetscPrintf(PETSC_COMM_WORLD," -- Reading fixations\n"); 
      // Read fixations
      // dofmap->fixa = da_create(sizeof(fixa_entry));
      // fixa_entry fe;
      nfixa=0;
      while (1) {
	fstack->get_line(line);
	astr_copy_s(linecopy, line);
	if (strstr(line,"__END_FIXA__")) break;
	nfixa++;
	int nr=sscanf(line,"%d %d %lf",&node,&kdof,&dval);
	if (nr !=3) {
	  PetscPrintf(PETSC_COMM_WORLD,
		      "Error reading fixations, for fixation %d\n",nfixa);
	  CHKERRQ(1);
	}

	dofmap->get_row(node,kdof,row);
	if (row.size()!=1) {
	  PetscPrintf(PETSC_COMM_WORLD,
		      "In line: %s\nFixation %d, imposed on an invalid node/field combination.\n",
		      astr_chars(linecopy),nfixa);
	}
	  
	int keq = row.begin()->first;
	assert(row.begin()->second == 1.);
	
	dofmap->fixed.push_back(fixation_entry(dval));
	fixed_dofs[keq]=dofmap->fixed.size()-1;

      }
      PetscPrintf(PETSC_COMM_WORLD,"Total fixations: %d\n",nfixa);

    } else if (!strcmp(token,"fixa_amplitude")) {

      PetscPrintf(PETSC_COMM_WORLD," -- Reading fixa_amplitude section\n"); 
      // next token is the identifier of the amplitude function
      token = strtok(NULL,bsp);
      char *label = new char[strlen(token)+1];
      strcpy(label,token);

      // create a new Amplitude and insert pointer in the list
      Amplitude *amp = new Amplitude(label);

      // amplitude_list:= stores a list of all the amplitudes defined 
      amplitude_list.push_back(amp);
      amp->read_hash_table(fstack);

      nfixa=0;
      while (1) {
	fstack->get_line(line);
	astr_copy_s(linecopy, line);
	if (strstr(line,"__END_FIXA__")) break;
	nfixa++;
	int nr=sscanf(line,"%d %d %lf",&node,&kdof,&dval);
	if (nr !=3) {
	  PetscPrintf(PETSC_COMM_WORLD,
		      "Error reading fixations, for fixation %d\n",nfixa);
	  CHKERRQ(1);
	}

	dofmap->get_row(node,kdof,row);
	if (row.size()!=1) {
	  PetscPrintf(PETSC_COMM_WORLD,
		      "In line: %s\nFixation %d, imposed on an invalid node/field combination.\n",
		      astr_chars(linecopy),nfixa);
	}
	  
	int keq = row.begin()->first;
	assert(row.begin()->second == 1.);
	
	dofmap->fixed.push_back(fixation_entry(dval,amp));
	fixed_dofs[keq]=dofmap->fixed.size()-1;

      }
      PetscPrintf(PETSC_COMM_WORLD,
		  "Total fixations with temporal amplitude: %d\n",nfixa);

    } else if (!strcmp(token,"constraint")) {

      PetscPrintf(PETSC_COMM_WORLD," -- Reading constraint section\n"); 
#define ERRLINE \
      if (ierr) { \
         PetscPrintf(PETSC_COMM_WORLD,"Error reading line \"%s\"",linecopy); \
         CHKERRQ(1); \
      }

      // Read constraints
      int node1,kdof1,nconstr=0;
      Constraint constraint;
      while (1) {
	
	fstack->get_line(line);
	if (strstr(line,"__END_CONSTRAINT__")) break;
	nconstr++;
	constraint.empty();
	double coef;
	int node,field;
	rflag=0;

	while (1) { 
	  ierr = readval(rflag,line,coef); if(ierr) break;
	  ierr = readval(rflag,line,node); ERRLINE;
	  ierr = readval(rflag,line,field); ERRLINE;
	  assert(node<=nnod);
	  assert(field<=ndof);
	  constraint.add_entry(node,field,coef);
	}
	dofmap->set_constraint(constraint);
      }
      PetscPrintf(PETSC_COMM_WORLD,"Total constraints: %d\n",nconstr);
#undef ERRLINE

    } else {
      PetscPrintf(PETSC_COMM_WORLD,"Bad section name in data file.\n"
		  "line: \"%s\"",line);
      CHKERRQ(1);
    }

  }

  // Read processor info
  const char *proc_weights;
  float *tpwgts;
  // Read data
  FileStack *weights_file;

  // tpwgts:= array containing the weight (speed) of the processor. To
  // be passed to METIS for a balanced partitioning.
  tpwgts = new float[size];
  dofmap->tpwgts = tpwgts; // add to dofmap
  mesh->global_options->get_entry("proc_weights",proc_weights);
  PetscPrintf(PETSC_COMM_WORLD,"size: %d\n",size);
  if (size==1 || proc_weights == NULL) {

    // If there is not  a processor weight table, then 
    for (int proc=0; proc<size; proc++) {
      tpwgts[proc] = 1./float(size);
    }

  } else {

    weights_file = new FileStack(proc_weights);
    float sumw=0.;
    for (int proc=0; proc<size; proc++) {
      weights_file->get_line(line);
      sscanf(line,"%f",&tpwgts[proc]);
      sumw += tpwgts[proc];
    }
    PetscPrintf(PETSC_COMM_WORLD,"total weight: %f\n",sumw);
    for (int proc=0; proc<size; proc++) {
      tpwgts[proc] /= sumw;
      PetscPrintf(PETSC_COMM_WORLD," proc: %d, w: %f\n",proc,tpwgts[proc]);
    }
    weights_file->close();
    
  }
    
  // nelemsets:= total number of elemesets in the mesh
  int nelemsets=da_length(mesh->elemsetlist);

  // Mesh partitioning
  int *npart = new int[nnod];
  if (size==1) {
    for (k=0; k<nnod; k++) {
      npart[k]=1;
    }
  }

  // is_fat:= flags if the elemset is fat or not. i.e. if it is an
  // elemset to be partitioned or not. Normally, we declare the
  // elemsets with large number of elements as fat.
  //
  // nelemfat:= total number of elements in fat elemsets
  // is_any_fat:= At least one elemset has to be fat
  // n2eptr,node2elem:= lists elements connected to a given element
  int nelemfat=0,is_fat,is_any_fat=0;
  int *n2eptr = new int[nnod+1];
  int *node2elem;
  int *nelemsetptr = new int[nelemsets+1];
  nelemsetptr[0]=0;

  for (node=0; node<=nnod; node++) 
    n2eptr[node]=0;

  // compute total number of elements in fat elemsets, nelemfat
  // also n2eptr[node] gets the total number of elements connected to
  // node and n2eptr[0]=0

  // numfat:= number of fat elemsets
  numfat=0;
  for (int ielset=0; ielset<nelemsets; ielset++) {
    elemset  = *(Elemset **)da_ref(mesh->elemsetlist,ielset);
    is_fat=1;
    ierr = get_int(elemset->thash,"is_fat",&is_fat,1); CHKERRQ(ierr);
    elemset->isfat=is_fat;
    nelemsetptr[ielset+1]=nelemsetptr[ielset];
    if (!is_fat) continue;
    numfat++;
    nelem = elemset->nelem;
    icone = elemset->icone;
    nel = elemset->nel;
    is_any_fat=1;
    nelemfat += nelem;
    nelemsetptr[ielset+1]=nelemsetptr[ielset]+nelem;
    for (int iel=0; iel<nelem; iel++) {
      for (int iloc=0; iloc<nel; iloc++) {
	// printf("adding to node %d\n",ICONE(iel,iloc));
	n2eptr[ICONE(iel,iloc)]++;
      }
    }
  }

  // raise error if not defined at least one fat elemset
  PFEMERRCQ(size>1 && !is_any_fat,
	    "At least one elemset has to be set as fat (is_fat property)\n"
	    "if numbre of processors is np>1\n");

  // cumulated sum 
  for (node=0; node<nnod; node++) 
    n2eptr[node+1] += n2eptr[node];

  // n2esize:= total number of entries in the node2elem table
  int n2esize = n2eptr[nnod];

  // initialize node2elem
  node2elem = new int[n2esize];
  for (int j=0; j<n2esize; j++) node2elem[j]=-1;

  // define node2elem pointer
  for (int ielset=0; ielset<nelemsets; ielset++) {
    elemset  = *(Elemset **)da_ref(mesh->elemsetlist,ielset);
    if (!elemset->isfat) continue;
    nelem = elemset->nelem;
    icone = elemset->icone;
    nel = elemset->nel;
    for (int iel=0; iel<nelem; iel++) {
      // iel:= global element number cumulated through the elemsetlist,
      // (only on fat elemsets)
      int ielg = nelemsetptr[ielset]+iel;
      for (int iloc=0; iloc<nel; iloc++) {
	// search for the next free position on the array
	node=ICONE(iel,iloc);
	int jpos;
	for (jpos=n2eptr[node-1]; jpos<n2eptr[node]; jpos++) 
	  if (node2elem[jpos]==-1) break;
	assert(jpos!=n2eptr[node]);
	node2elem[jpos]=ielg;
      }
    }
  }

#if 0
  for (node=0; node<nnod; node++) {
    printf("node %d: ",node+1);
    for (int jpos=n2eptr[node]; jpos<n2eptr[node+1]; jpos++) 
      printf(" %d",node2elem[jpos]);
    printf("\n");
  }
#endif

  // wait_from_console("before partitioning"); 
  PetscPrintf(PETSC_COMM_WORLD,"Starts partitioning.\n"); 

#define USE_METIS_PARTITION
#ifdef USE_METIS_PARTITION
  // Create adjacency table for aprtitioning with Metis. In the adjacency
  // graph the nodes are the elements of the FEM mesh. Two nodes of
  // the graph (elements of the mesh) have are linked if they share a
  // node. 

  // adjncy:= xadj:= graph desdcribed in CSR format (as defined in
  // Metis documentation)

  int *vpart = new int[nelemfat];
  { // This artificial block is in order to adjncy be deleted.
    vector<int> adjncy;
    int *xadj = new int[nelemfat+1];

    // mark:= auxiliary vector that flags if an element has been marked
    // already as a linked node in the graph
    int *mark = new int[nelemfat];
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
	      adjncy.push_back(ielgjj);
	      mark[ielgjj]=ielgj;
	      xadj[ielgj+1]++;
	      // printf("adding %d %d\n",ielgj,ielgjj);
	    }
	  }
	}
      }
    }

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
      METIS_WPartGraphKway(&nelemfat,xadj,adjncy.begin(),NULL, 
			   NULL,&wgtflag,&numflag,&size, 
			   tpwgts,&options,&edgecut,vpart);
    } else {
      for (int jj=0; jj<nelemfat; jj++) 
	vpart[jj]=0;
    }

    // wait_from_console("antes de borrar cosas para part.");  
    delete[] xadj;
    VOID_IT(adjncy);
    // adjncy.~vector();
    delete[] mark;
    // delete[] tpwgts; // don't delete, keep it in dofmap
    // wait_from_console("despues de borrar cosas para part.");  
  }
  // wait_from_console("despues de borrar adjncy");  

#else //USE_METIS_PARTITION
  int *vpart = new int[nelemfat];
  if (size>1) {
    for (int jj=0; jj<nelemfat; jj++) {
      vpart[jj]= jj/(nelemfat/size+1);
    }
  } else {
    for (int jj=0; jj<nelemfat; jj++) 
      vpart[jj]=0;
  }
#endif

  PetscPrintf(PETSC_COMM_WORLD,"Ends partitioning.\n");

  // nelem_part:= nelem_part[proc] is the number of elements in
  // processor proc. 
  int *nelem_part = new int[size];
  for (int proc=0; proc<size; proc++) nelem_part[proc]=0;
  for (int jj=0; jj<nelemfat; jj++) {
    nelem_part[vpart[jj]]++;
    // printf("elem %d in proc %d\n",jj+1,vpart[jj]);
  }
  if (size>1 && myrank==0) {
    for (int proc=0; proc<size; proc++) 
      printf("%d elements in processor %d\n",nelem_part[proc],proc);
  }

  // Partition nodes. Assign to each node the partition of the first
  // element in tne node2elem list. If there is no element, then
  // assign processor 0 and give a warning.

  // In the future this may be done in a better way. Following links
  // in the other elemensets, even if these elemensets have not been
  // taken into account in the partitioning. 

  // node_not_connected_to_fat:= flags whether there is a node not
  // connected to any fat elemset or not. 
  int node_not_connected_to_fat=0;
  
  // wait_from_console("trace 0");  

  int print_nodal_partitioning=0;
  ierr = get_int(mesh->global_options,
		 "print_nodal_partitioning",
		 &print_nodal_partitioning,1); CHKERRQ(ierr);
  if (print_nodal_partitioning)
    PetscPrintf(PETSC_COMM_WORLD,"\nNodal partitioning (node/processor): \n");
  for (node=0; node<nnod; node++) {
    if (n2eptr[node]==n2eptr[node+1]) {
      node_not_connected_to_fat=1;
      npart[node]=1;
    } else {
      npart[node] = vpart[node2elem[n2eptr[node]]]+1;
    }
    if (print_nodal_partitioning)
      PetscPrintf(PETSC_COMM_WORLD,
		  "%d   %d\n",node+1,npart[node]);
  }
  if (print_nodal_partitioning)
    PetscPrintf(PETSC_COMM_WORLD,"End nodal partitioning table\n\n");

  if (node_not_connected_to_fat)
    PetscPrintf(PETSC_COMM_WORLD,"warning! there is at least one"
		" node not linked to any fat elemset. \n"
		"This induces artificial numbering\n");
  
  // define the eparts of each fat elemset as the corresponding part
  // of vpart. 
  GETOPTDEF(int,debug_element_partitioning,0);
  for (int ielset=0; ielset<nelemsets; ielset++) {
    elemset  = *(Elemset **)da_ref(mesh->elemsetlist,ielset);
    if (!elemset->isfat) continue;
    int *epart= new int[elemset->nelem];
    elemset->epart = epart;
    for (int iel=0; iel<elemset->nelem; iel++) {
      int ielg=nelemsetptr[ielset]+iel;
      epart[iel] = vpart[ielg]+1;
    }
    if (debug_element_partitioning) {
      PetscPrintf(PETSC_COMM_WORLD,
		  "Element \"%s\", type \"%s\", nelem %d\n."
		  "Partitioning table: \n");
      for (int kk=0; kk<nelem; kk++) 
	PetscPrintf(PETSC_COMM_WORLD,"%d -> %d\n",kk,epart[kk]);
    }
  }
  
  if (debug_element_partitioning==2) {
    PetscFinalize();
    exit(0);
  }

  // wait_from_console("trace 1");  

  delete[] n2eptr;
  delete[] node2elem;
  delete[] nelemsetptr;
  delete[] vpart;
  delete[] nelem_part;

  // wait_from_console("trace 2");  

  // Number degrees of freedom
  int proc,*neqproc, *startproc;
  neqproc= new int[size+1];
  startproc= new int[size+1];
  int *perm = new int [ndof*nnod];
  for (k=0; k<ndof*nnod; k++) {
    perm[k]=0;
  }

  // wait_from_console("trace 3");  

  // First number all the edof's in processor 0, then on 1, etc...
  // perm contains the permutation. 

  // The criterion to assign a column index (dof or fixed) to a
  // processor is as follows: take the corresponding column, take the
  // first non-null entry, then assign to the processor that the
  // corresponding node belongs. 

  int field;
  map<int,int>::iterator end_fix;
  end_fix = fixed_dofs.end();
  
  // First, put perm = to minus the corresponding processor. 
  // The fixed dof's are set to `-(size+1)'  (number of processors + 1)
  for (k=1; k<=nnod*ndof; k++) {
    dofmap->id->get_col(k,col);
    if (col.size()==0) continue;
    if (fixed_dofs.find(k) != end_fix) {
      perm[k-1] = -(size+1);
      continue;
    }
    edof = col.begin()->first;
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
  PetscPrintf(PETSC_COMM_WORLD,
	      "Total number of degrees of freedom neq:     %d\n"
	      "Total number of independent fixations neqf: %d\n",
	      neq,dofmap->neqf);

#if 0
  // This is the old version related to the remap_cols() bug. A leak
  // of memory due to the memory management of the STL map routines. 
  dofmap->id->remap_cols(perm);
#else
  idmap *idnew = new idmap(nnod*ndof,NULL_MAP);
  dofmap->id->remap_cols(perm,*idnew);
  delete dofmap->id;
  dofmap->id = idnew;
#endif
 
  TGETOPTDEF(mesh->global_options,int,check_dofmap_id,0);
  if (check_dofmap_id && myrank==0) {
    dofmap->id->check();
    dofmap->id->print_statistics();
  }
  TGETOPTDEF(mesh->global_options,int,print_dofmap_id,0);
  if (print_dofmap_id && myrank==0) {
    dofmap->id->print("dofmap->id: ");
  }

  map<int,int>::iterator jj;

#if 0  // debug:=
  {
    map<int,int>::iterator jjjj;
    for (jjjj=fixed_dofs.begin(); jjjj!=end_fix; jjjj++) 
      printf("%d -> fix %d\n",jjjj->first,jjjj->second);

    for (int kkk=0; kkk<nnod*ndof; kkk++) {
      printf("%d -> perm[%d]\n",kkk,perm[kkk]);
    }
  
  }
#endif

  int nfixed = dofmap->fixed.size();
  vector<fixation_entry> fixed_remapped(nfixed);
  //  fixed_remapped.reserve(nfixed);
  for (jj=fixed_dofs.begin(); jj!=fixed_dofs.end(); jj++) {
    int newjeq = perm[jj->first -1];
    if (newjeq==0) continue;
    int indx= newjeq - neq - 1;
    assert(0<=indx && indx<nfixed); //This should be in this range. Otherwise it falls
				    // off the `fixed_remapped' vector.
    fixed_remapped[indx] =
      (dofmap->fixed)[jj->second];
  }

  // wait_from_console("trace 6");  

  // swap fixed with fixed_remapped (reordered)
  dofmap->fixed.swap(fixed_remapped);
  VOID_IT(fixed_remapped);
  VOID_IT(fixed_dofs);

#if 0
  //debug:= Verificar como quedaron las fijaciones
  for (int j=dofmap->neq+1; j<=dofmap->neqtot; j++) {
    dofmap->id->get_col(j,col);
    assert(col.size()==1);
    edof = col.begin()->first;
    assert(col.begin()->second == 1.);
    dofmap->nodf(edof,node,field);
    printf("node %d, field %d , edof %d, j %d, val %f\n",
	   node, field , edof,j,fixed_remapped[j-neq-1]);
  }
#endif
  delete[] perm;

  // Dof's connected to elements in this processor
  int jel,locdof,keq;
  dof_here = new int [neq];

  for (keq=0; keq<neq; keq++) dof_here[keq]=0;

  int dof1,dof2; // interval of dof's that live in the processor
  dof1=startproc[myrank]+1;
  dof2=dof1+neqproc[myrank]-1;
  dofmap->dof1 = dof1;
  dofmap->dof2 = dof2;
  set<int> ghost_dof_set;
    
  // wait_from_console("trace 7");  

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
	if(elemset->epart[iele]!=myrank+1) continue;
	elemset->nelem_here++;
	for (jel=0; jel<nel; jel++) {
	  node = ICONE(iele,jel);
	  for (kdof=1; kdof<=ndof; kdof++) {
	    dofmap->get_row(node,kdof,row);
	    for (kndx=row.begin(); kndx!=row.end(); kndx++) {
	      int dof = kndx->first;
	      if (dof <= neq)
		dof_here[dof-1]=1;
	      if (dof <= neq && !(dof1 <= dof && dof <= dof2))
		ghost_dof_set.insert(dof-1);
	    }
	  }
	}
      }
    } else {
      elemset->epart = new int[nelem];
      for (iele=0; iele<nelem; iele++) {
	// Decide in which processor will be computed this element
	node = ICONE(iele,0);
	proc = npart[node-1];
	if (proc<1 || proc>size) {
	  PetscPrintf(PETSC_COMM_WORLD,
		      "node %d attached to processor %d out of range",node,proc);
	  CHKERRA(1);
	}
	elemset->epart[iele]=proc;

	if(elemset->epart[iele]!=myrank+1) continue;
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
    // ghost_elems:= These are elements that have related dof's on the
    // processor, but they don't live on the processor.
    // Loops for defining profiles of matrices must loop over them
    // also. 

    // wait_from_console("trace 8");  

    elemset->ghost_elems = da_create(sizeof(int));

    for (iele=0; iele<nelem; iele++) {
      if(elemset->epart[iele]==myrank+1) continue;
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
      elemset->local_store = new (void *)[elemset->nelem];
      for (int j=0; j<elemset->nelem_here; j++) {
	elemset->local_store[j]=NULL;
      }
    }

    da_sort (elemset->ghost_elems,int_cmp,NULL);
    int nghostel = da_length(elemset->ghost_elems);

    if (size>1) {
      PetscPrintf(PETSC_COMM_WORLD,"For elemset \"%s\"\n",elemset->type);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			      "On processor [%d], %d local elements,"
			      " %d ghost elements.\n",
			      myrank,elemset->nelem_here,nghostel);
      PetscSynchronizedFlush(PETSC_COMM_WORLD);
    }

  }
  
  // wait_from_console("trace 10");  

  // Convert ghost_dof_set en dofmap
  set<int>::iterator it;
  dofmap->ghost_dofs = new vector<int>;
  for (it=ghost_dof_set.begin(); it!=ghost_dof_set.end(); it++) 
    dofmap->ghost_dofs->push_back(*it);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  VOID_IT(ghost_dof_set);
  int nghost_dofs = dofmap->ghost_dofs->size();

  // This dof_here stuff may be done with STL sets
  int ndofhere = 0;
  for (keq=0; keq<neq; keq++) {
    ndofhere += dof_here[keq];
  }

  int *dof_here_list;
  dof_here_list = new int[ndofhere];

  k=0;
  for (keq=0; keq<neq; keq++) {
    if (dof_here[keq]) {
      dof_here_list[k++] = keq;
    }
  }
  // fixme:= Al final ya practicamente saque todo lo de dof_here y
  // dof_here_list, asi que esto de arriba y bajo habria que sacarlo. 
  delete[] dof_here;
  delete[] dof_here_list;

  // wait_from_console("trace 11");  

  // Defines certain quantities in dofmap
  dofmap->neq= neq;
  dofmap->startproc = startproc;
  dofmap->neqproc = neqproc;
  dofmap->size = size;
  dofmap->npart = npart;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  //                        DEFINE SCATTER
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

  // wait_from_console("trace 12");  

  PetscPrintf(PETSC_COMM_WORLD,"Defining scatters...\n");

  // Create MPI vectors
  //  Vec x,xseq;
  Vec x,ghost_vec,xseq;
  IS is_ghost_glob,is_ghost_loc,is_print;

  dofmap->ghost_scatter = new VecScatter;
  dofmap->scatter_print = new VecScatter;

  dofmap->create_MPI_vector(x);
  dofmap->create_MPI_ghost_vector(ghost_vec);

  int nghost_tot;
  ierr = VecGetSize(ghost_vec,&nghost_tot); CHKERRQ(ierr);

  ierr = ISCreateGeneral(PETSC_COMM_WORLD,nghost_dofs,
			 dofmap->ghost_dofs->begin(),
			 &is_ghost_glob);  CHKERRQ(ierr); 
  ierr = ISCreateStride(PETSC_COMM_WORLD,nghost_dofs,0,1,&is_ghost_loc);

  ierr = VecScatterCreate(x,is_ghost_glob,ghost_vec,is_ghost_loc,
			  dofmap->ghost_scatter); CHKERRQ(ierr); 

  ierr = ISDestroy(is_ghost_glob); CHKERRQ(ierr); 
  ierr = ISDestroy(is_ghost_loc); CHKERRQ(ierr); 

  int neql = (myrank==0 ? dofmap->neq : 0);
  ierr = VecCreateSeq(PETSC_COMM_SELF,neql,&xseq);  CHKERRQ(ierr);
  
  ierr = ISCreateStride(PETSC_COMM_WORLD,neql,0,1,&is_print);
  ierr = VecScatterCreate(x,is_print,xseq,is_print,
			  dofmap->scatter_print); CHKERRQ(ierr); 
  
  ierr = ISDestroy(is_print); CHKERRQ(ierr); 

#if 0
  for (int jj=dof1; jj<=dof2; jj++) {
    VecSetValue(x,jj-1,double(jj),INSERT_VALUES);
  }
  PetscPrintf(PETSC_COMM_WORLD,"vector x en read_mesh\n");
  ierr = VecView(x,VIEWER_STDOUT_WORLD); CHKERRA(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"================\n");
  ierr = VecScatterBegin(x,ghost_vec,INSERT_VALUES,
			 SCATTER_FORWARD,dofmap->ghost_scatter); CHKERRA(ierr); 
  ierr = VecScatterEnd(x,ghost_vec,INSERT_VALUES,
		       SCATTER_FORWARD,dofmap->ghost_scatter); CHKERRA(ierr); 
  
  double *array;
  ierr = VecGetArray(ghost_vec,&array); CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Local ghost values on [%d]\n",myrank);
  for (int jj=0; jj<nghost_dofs; jj++ ) 
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"local %d, dof %d  -> %g\n",
			    jj,(*dofmap->ghost_dofs)[jj],array[jj]);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
#endif

  ierr = VecDestroy(x);
  ierr = VecDestroy(ghost_vec);
  ierr = VecDestroy(xseq);

#if 0
  ierr = VecGhostGetLocalForm(gx,&lx); CHKERRQ(ierr);

  for (int jj=dof1; jj<=dof2; jj++) {
    VecSetValue(gx,jj-1,double(jj),INSERT_VALUES);
  }
  ierr = VecAssemblyBegin(gx); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(gx); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(gx,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(gx,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  double *array;
  ierr = VecGetArray(lx,&array); CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"On processor [%d], local values\n",myrank);
  for (int jj=0; jj<neqproc[myrank]; jj++ ) {
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"local %d, global %d -> %g\n",
			    jj,startproc[myrank]+jj,array[jj]);
  }

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Local ghost values\n");
  for (int jj=0; jj<nghost_dofs; jj++ ) {
    int local = neqproc[myrank]+jj;
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"local %d, global %d -> %g\n",
			    local,(*(dofmap->ghost_dofs))[jj],array[local]);
  }

  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  ierr = VecRestoreArray(lx,&array);CHKERRQ(ierr);
  ierr = VecGhostRestoreLocalForm(gx,&lx);CHKERRQ(ierr); 
  PetscFinalize();
  exit(0);
#endif 

#if 0
  // para debuggear
  ierr = VecScatterBegin(x,xseq,INSERT_VALUES,SCATTER_FORWARD,scatter); CHKERRQ(ierr); 
  ierr = VecScatterEnd(x,xseq,INSERT_VALUES,SCATTER_FORWARD,scatter); CHKERRQ(ierr); 
  PetscPrintf(PETSC_COMM_WORLD,"Despues del scatter\n");


  // debug
  ierr = VecGetArray(xseq,&sstate); CHKERRQ(ierr); 
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"On processor %d ------------------\n",
			  myrank+1);
  for (int k=0; k<neq; k++) {
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d -> %f\n",k+1,sstate[k]);
  }
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  ierr = VecRestoreArray(xseq,&sstate); CHKERRQ(ierr); 

  PetscFinalize();
  exit(0);
#endif  

#if 0
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			  "On processor %d \n",myrank+1);
  for (k=0; k<ndofhere; k++) {
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			    "%d -> %d\n",k+1,dof_here_list[k]);
  }
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
#endif
  
  PetscPrintf(PETSC_COMM_WORLD,"Total number of dof's: %d\n",neq);

#if 0
  int nloads=0;
  fstack->get_line(line);
  if (strncmp("loads",line,5) ) {
    PetscPrintf(PETSC_COMM_WORLD,"Couldn't find <fixa> tag line\n");
    exit(1);
  }
  while (1) {
    fstack->get_line(line);
    if (strstr(line,"__END_LOADS__")) break;
    nloads++;
    int nr=sscanf(line,"%d %d %lf",&node,&kdof,&dval);
    if (nr !=3) {
      PetscPrintf(PETSC_COMM_WORLD,"Error reading LOADS, for loads %d\n",nfixa);
    }
    //LOAD(node-1,kdof)=dval;
    CHKERRQ(1);
  }
  PetscPrintf(PETSC_COMM_WORLD,"Total loads: %d\n",nloads);
#endif

  //fclose(fid);
  astr_destroy(linecopy);
 
  return 0;

}
