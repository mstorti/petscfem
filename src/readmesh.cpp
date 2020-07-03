//__INSERT_LICENSE__
//$Id merge-with-petsc-233-55-g52bd457 Fri Oct 26 13:57:07 2007 -0300$
#ifndef _GNU_SOURCE 
#define _GNU_SOURCE 
#endif
#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/idmap.h>
#include <src/idmap.h>
#include <src/generror.h>
#include <src/syncbuff.h>
#include <src/debug.h>

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
#include <queue>
#include <src/getprop.h>
#include <src/pfobject.h>
#include <src/dvecpar.h>
#include <src/mshpart.h>
#include <src/h5utils.h>

// #define TRACE_READMESH
#ifdef TRACE_READMESH
#undef TRACE
#define TRACE(n)				\
   GLOBAL_DEBUG->trace("readmesh trace " #n);
#else
#undef TRACE
#define TRACE(n)
#endif

using namespace std;
Mesh *GLOBAL_MESH;
#define IDENT(j,k) VEC2(ident,j,k,ndof)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
static int ishdf5(const char *filedset) {
  size_t len = strlen(filedset);
  for (size_t j=0; j<len; j++)
    if (filedset[j]==':') return 1;
  return 0;
}

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
#undef ICONE
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define ELEMIPROPS(j,k) VEC2(elemiprops,j,k,neliprops)
#define NODEDATA(j,k) VEC2(mesh->nodedata->nodedata,j,k,nu)

#define RM_ASSERT(cond,mess)					\
PETSCFEM_ASSERT(cond,mess "%s:%d: at (or after) line: \"%s\"",	\
fstack->file_name(),fstack->line_number(),fstack->line_read());

// fixme:= aca no se porque tuve que pasar neq por referenciar
// porque sino no pasaba correctamente el valor. 
#undef __FUNC__
#define __FUNC__ "read_mesh"
int read_mesh(Mesh *& mesh,char *fcase,Dofmap *& dofmap,
	      int & neq,int size,int myrank) {
  vector<Amplitude *> amplitude_list;

  char *p1,*p2, *token, *type;
  const char *bsp=" \t\n";
  vector<string> tokens;
  char *line;
  const char *cline;
  Autostr *linecopy = astr_create();
  char *key,*val;
  int ndim,nu,ndof,nnod=0,ierr,ierro,numfat,node,jdof,kdof,edof;
  int pos, nelem, nel, nelprops, neliprops, nread, elemsetnum=0,
	iele,k,nfixa, *ident,rflag;
  double *dptr,dval; 

  TextHashTable *thash=NULL;
  Nodedata *nodedata=NULL;

  row_t row,col,row0;
  row_t::iterator kndx;

  Elemset *elemset=NULL;
  int *icone=NULL,etype,*dof_here=NULL;
  mesh = new Mesh;
  GLOBAL_MESH = mesh;
  mesh->nodedata = new Nodedata;
  mesh->elemsetlist = da_create(sizeof(Elemset *));

  // Read data
  FileStack *fstack;
  fstack = new FileStack(fcase);
  fstack->quiet = !myrank;	// Only output on server
  PETSCFEM_ASSERT(fstack->ok(),"read_mesh: "
		  "Couldn't open file \"%s\"\n",fcase);
  if (myrank==0) fstack->set_echo_stream(stdout);

#ifdef TRACE_READMESH
    GLOBAL_DEBUG->activate("print");
    Debug::init();
#endif

  // Initialize number of eqs.
  neq=0;
  dofmap = new Dofmap;
      
  while (!fstack->get_line(line)) {
      
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
      
      TRACE(-5.0);
      PetscPrintf(PETSCFEM_COMM_WORLD," -- Reading nodes:\n");
      tokens.clear();
      tokenize(fstack->line_read(),tokens);
      assert(tokens.size()==4);
      ierr = string2int(tokens[1],ndim);
      assert(!ierr);
      assert(ndim>0);
      ierr = string2int(tokens[2],nu);
      assert(!ierr);
      assert(nu>=ndim);
      ierr = string2int(tokens[3],ndof);
      assert(!ierr);
      assert(ndof>=0);

      PetscPrintf(PETSCFEM_COMM_WORLD, 
		  "Dimension: %d, Size of nodedata vector: %d\n",ndim,nu);
      mesh->nodedata->nu = nu;
      mesh->nodedata->ndim = ndim;

      dofmap->ndof = ndof;
      node = 0;
      
      double *row = new double[nu];
      darray *xnod;
      xnod = da_create(nu*sizeof(double));
      while (!fstack->get_line(line)) {
	astr_copy_s(linecopy, line);
	if (strstr("__END_NODES__",line)) break;
	node++;
	for (int kk=0; kk<nu; kk++) {
	  token = strtok((kk==0 ? line : NULL),bsp);
	  
	  RM_ASSERT(token!=NULL,"Error reading coordinate[1]s\n");
	  int nread = sscanf(token,"%lf",row+kk);
	  RM_ASSERT(nread==1,"Error reading coordinates[2]\n");
	}
	int indx = da_append (xnod,row);
	if (indx<0) PFEMERRQ("Insufficient memory reading nodes");
      }
      RM_ASSERT(fstack->last_error()==FileStack::read_ok,
		"Error reading coordinates[3]\n");

      nnod=node;
      dofmap->nnod = nnod;
      mesh->nodedata->nnod = nnod;
      PetscPrintf(PETSCFEM_COMM_WORLD,"Read %d nodes\n",nnod);
      
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

    } else if (!strcmp(token,"nodedata")) {

      int nread;

      PetscPrintf(PETSCFEM_COMM_WORLD," -- Reading nodes:\n");

#define RM_READ_INT(var)                                                \
      token = strtok(NULL,bsp);                                         \
      PETSCFEM_ASSERT(token,#var " is required! line read: \"%s\"",         \
      fstack->line_read());                                             \
      nread = sscanf(token,"%d",&var);                                  \
      PETSCFEM_ASSERT(nread==1,"Could'nt read " #var "  in line: \"%s\", token \"%s\"", \
                      fstack->line_read(),token); 

      RM_READ_INT(ndim);
      RM_READ_INT(nu);
      RM_READ_INT(ndof);
 
      PetscPrintf(PETSCFEM_COMM_WORLD, 
		  "Dimension: %d, Size of nodedata vector: %d\n",ndim,nu);
      read_hash_table(fstack,mesh->nodedata->options);
      mesh->nodedata->nu = nu;
      mesh->nodedata->ndim = ndim;

      dofmap->ndof = ndof;
      node = 0;
      double *row = new double[nu];
      darray *xnod=NULL;
      const char *data = NULL;
      mesh->nodedata->options->get_entry("data",data);
      assert(data);

      thash = mesh->nodedata->options;
      TGETOPTDEF(thash,int,use_hdf5,0);
      TGETOPTDEF(thash,int,use_hdf5_auto,0);
      int hdf5 = use_hdf5 || (use_hdf5_auto && ishdf5(data));
      if (hdf5) {
        TGETOPTDEF_S(thash,string,data,NONE);
        PETSCFEM_ASSERT0(data!="NONE","data entry is required if use_hdf5");  
 
        dvector<double> dvxnod;
        if (!myrank) {
          // Read the vector in the master
          h5_dvector_read(data.c_str(),dvxnod);
          dvxnod.defrag();
          printf("read %d doubles\n",dvxnod.size());
          nnod = dvxnod.size();
          PETSCFEM_ASSERT(nnod%nu==0,
                          "bad size HDF5 xnod size %d",nnod);
        }
        ierr = MPI_Bcast(&nnod,1,MPI_INT,0,PETSCFEM_COMM_WORLD);
        CHKERRQ(ierr);
        nnod /= nu;
        mesh->nodedata->nodedata = new double[nnod*nu];
        dofmap->nnod = nnod;
        mesh->nodedata->nnod = nnod;
        if (!myrank) {
          for (node=0; node<nnod; node++) {
            for (int kk=0; kk<nu; kk++) 
              NODEDATA(node,kk) = dvxnod.e(node,kk);
          }
        }
        dvxnod.clear();
        ierr = MPI_Bcast (mesh->nodedata->nodedata,nnod*nu,
                          MPI_DOUBLE,0,PETSCFEM_COMM_WORLD);
      } else {
        AutoString datas;
        datas.set(data).deblank();
        data = datas.str();
        FileStack *fstack_nodes_data=NULL;
        if (!myrank) {
          fstack_nodes_data = new FileStack(data);
          ierro = !fstack_nodes_data->ok();
        }
        CHECK_PAR_ERR(ierro,"Error reading nodes");
        if (!myrank) try {
            xnod = da_create(nu*sizeof(double));
            while (!fstack_nodes_data->get_line(line)) {
              astr_copy_s(linecopy, line);
              node++;
              for (int kk=0; kk<nu; kk++) {
                token = strtok((kk==0 ? line : NULL),bsp);
                PETSCFEM_ASSERT_GE1(token!=NULL,
                                    "Error reading coordinates in line[4]:\n\"%s\"\n"
                                    "Not enough values in line!!\n",astr_chars(linecopy));
                int nread = sscanf(token,"%lf",row+kk);
                PETSCFEM_ASSERT_GE1(nread == 1,
                                    "Error reading coordinates in line[5]:\n\"%s\"",line);
              }
              int indx = da_append (xnod,row);
              PETSCFEM_ASSERT_GE0(indx>=0,"Insufficient memory reading nodes");
            } 
            ierro = fstack_nodes_data->last_error()!=FileStack::eof;
            PETSCFEM_ASSERT_GE1(!ierro,"Couldn't process correctly node data file %s\n",
                                fstack_nodes_data->file_name());
            nnod=node;
            fstack_nodes_data->close();
            delete fstack_nodes_data;
          } 
        CHECK_PAR_ERR_GE;

        ierr = MPI_Bcast (&nnod,1,MPI_INT,0,PETSCFEM_COMM_WORLD);
        dofmap->nnod = nnod;
        mesh->nodedata->nnod = nnod;
        PetscPrintf(PETSCFEM_COMM_WORLD,"Read %d nodes\n",nnod);
      
        delete[] row;
        mesh->nodedata->nodedata = new double[nnod*nu];
        if (!myrank) {
          for (node=0; node<nnod; node++) {
            row = (double *) da_ref(xnod,node);
            for (int kk=0; kk<nu; kk++) {
              NODEDATA(node,kk) = row[kk];
            }
          }
          da_destroy(xnod);
        }

        ierr = MPI_Bcast (mesh->nodedata->nodedata,nnod*nu,
                          MPI_DOUBLE,0,PETSCFEM_COMM_WORLD);
      }

      // calling dofmap constructor
      dofmap->id = new idmap(nnod*ndof,NULL_MAP);

    } else if (!strcmp(token,"table")) {

      token = strtok(NULL,bsp);
      PETSCFEM_ASSERT0(token,"Couldn't find name for table section");  
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
      token = strtok(NULL,bsp);
      PETSCFEM_ASSERT0(token,"Couldn't find elemset name.");  
      type = new char[strlen(token)+1];
      strcpy(type,token);
      PetscPrintf(PETSCFEM_COMM_WORLD,"\n -- Reading elemset type \"%s\"\n",type);
      token = strtok(NULL,bsp); 
      PETSCFEM_ASSERT0(token,"Couldn't find elemset `nel' parameter.");  
      int nread = sscanf(token,"%d",&nel);
      PETSCFEM_ASSERT0(nread==1,"Couldn't convert `nel' parameter to an integer.");  

      // Reads hash table for the elemset
      read_hash_table(fstack,thash);

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
	for (unsigned int j=0; j<prop_names.size(); j++) {
	  prop_name = prop_names[j].c_str();
	  phe = new props_hash_entry;
	  // fixme:= despues hay que poner width a la cantidad de
	  // escalares que contiene esta propiedad
	  phe->width=prop_lens[j];
	  phe->position = posit;
	  posit += phe->width;

	  // verify that property is not duplicated
	  // fixme:= cast to `(char *)' is for avoiding a warning with
	  // old compiler versions
	  pheold = (props_hash_entry *)g_hash_table_lookup(props,(char *)prop_name);
	  if (pheold) {
	    PetscPrintf(PETSCFEM_COMM_WORLD,
			"duplicated elem properties label \"%s\"\n",prop_name);
	    PetscFinalize();
	    exit(1);
	  }
	  key = new char[strlen(prop_name)+1];
	  strcpy(key,prop_name);
	  g_hash_table_insert (props,(void *)key,phe);
	}
	// g_hash_table_freeze(props);
      }
      nelprops = posit;

      // Read iprops line
      thash->get_entry("iprops",cline);
      neliprops=0; 
      // GHashTable *iprops = g_hash_table_new(&g_str_hash,&g_str_equal);
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
	    PetscPrintf(PETSCFEM_COMM_WORLD,
			"duplicated elem properties label \"%s\"\n",token);
	    PetscFinalize();
	    exit(1);
	  }
	  key = new char[strlen(token)+1];
	  strcpy(key,token);
	  g_hash_table_insert (props,(void *)key,phe);
	}
	// g_hash_table_freeze(props);
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

      int nelprops_add=-1, neliprops_add=-1;
      double *elemprops_add=NULL;
      int *elemiprops_add=NULL;
      double *elemprops = NULL;
      int *elemiprops = NULL;

      TGETOPTDEF_S(thash,string,data,NONE);
      TGETOPTDEF(thash,int,use_hdf5,0);
      TGETOPTDEF(thash,int,use_hdf5_auto,0);
      int hdf5 = use_hdf5 || (use_hdf5_auto && ishdf5(data.c_str()));
      if (hdf5) {
        // PETSCFEM_ASSERT(size==1,"Not implemented yet MPI size>1. Size=%d",size);
        PETSCFEM_ASSERT0(data!="NONE","data entry is required if use_hdf5");  
        // TGETOPTDEF_S(thash,string,dset,NONE);
        // PETSCFEM_ASSERT0(dset!="NONE","dset entry is required if use_hdf5");  
        dvector<int> dvicone;
        if (!myrank) {
          h5_dvector_read(data.c_str(),dvicone);
          dvicone.defrag();
          printf("read %d ints\n",dvicone.size());
          nelem = dvicone.size();
          PETSCFEM_ASSERT(nelem%nel==0,
                          "bad elemset HDF5 dataset size, "
                          "size %d, nel %d",nelem,nel);  
          nelem /= nel;
          dvicone.reshape(2,nelem,nel);
        }
        ierr = MPI_Bcast(&nelem,1,MPI_INT,0,PETSCFEM_COMM_WORLD);

        elemprops = new double[nelem*nelprops];
        elemiprops = new int[nelem*neliprops];
        icone = new int[nelem*nel];

        dvector_clone_parallel(dvicone);
        dvicone.defrag();
        for (iele=0; iele<nelem; iele++) {
          for (int kk=0; kk<nel; kk++) {
            node = dvicone.e(iele,kk);
            ICONE(iele,kk) = node;
            for (int kdof=1; kdof<=ndof; kdof++) {
              edof = dofmap->edof(node,kdof);
              dofmap->id->set_elem(edof,edof,1.);
            }
          }
        }
        dvicone.clear();
        ierr = MPI_Bcast(icone,nelem*nel,MPI_INT,0,PETSCFEM_COMM_WORLD);

        TGETOPTDEF(thash,int,additional_props,0);
        nelprops_add = additional_props;
        TGETOPTDEF(thash,int,additional_iprops,0);
        neliprops_add = additional_iprops;

        // Reading per element property  w/HDF5
        if (nelprops>0) {
          if (!myrank) {
            dvector<double> dvprops;
            TGETOPTDEF_S(thash,string,dataprops,NONE);
            PETSCFEM_ASSERT(dataprops!="NONE",
                            "dataprops entry is required if use_hdf5 and "
                            "nelprops>0. nelprops: %d",nelprops);  
            h5_dvector_read(dataprops.c_str(),dvprops);
            dvprops.defrag();
            // printf("dataprops read %d doubles\n",dvprops.size());
            int sz = dvprops.size();
            PETSCFEM_ASSERT(sz==nelem*nelprops,
                            "Bad dataprops size %d, nelem %d, nelprops %d",
                            sz,nelem,nelprops);
            dvprops.reshape(2,nelem,nelprops);
            // reading element properties
            for (int k=0; k<nelem; k++) 
              for (int jprop=0; jprop<nelprops; jprop++) 
                ELEMPROPS(k,jprop) = dvprops.e(k,jprop);
            dvprops.clear();
          }
          ierr = MPI_Bcast(elemprops,nelem*nelprops,MPI_DOUBLE,0,
                           PETSCFEM_COMM_WORLD);
        }

        // Check that other 'per element' properties are
        // void, because until it is implemented in HDF5
#define CHECKVAR(var) PETSCFEM_ASSERT(var==0,                      \
          "Not implemented yet for HDF5 connectivities. "       \
          "%s %d",#var,var);
        CHECKVAR(neliprops);
        CHECKVAR(nelprops_add);
        CHECKVAR(neliprops_add);
      } else {
        const char *data=NULL;
        thash->get_entry("data",data);
        if (!data) {
          while (1) {
            ierr = fstack->get_line(line);
            PETSCFEM_ASSERT0(ierr==0,"Couldn't find __END_ELEMSET__ tag");  
	  
            // reading element connectivities
            for (int jel=0; jel<nel; jel++) {
              token =  strtok(( jel==0 ? line : NULL),bsp);
              if (!token) fstack->print();
              PETSCFEM_ASSERT(token,
                              "Error reading element connectivities at\n"
                              "%s:%d: at (or after) line: \"%s\"",fstack->file_name(),
                              fstack->line_number(),
                              fstack->line_read());
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
                PetscPrintf(PETSCFEM_COMM_WORLD,
                            "fails to read per-element property %d,\n at line \"%s\"",
                            jprop+1,line);
                PFEMERRQ("\n");
              }
              sscanf(token,"%lf",proprow+jprop);
            }

            // reading integer element properties
            for (int jprop=0; jprop<neliprops; jprop++) {
              token = strtok(NULL,bsp);
              if (token==NULL) {
                PetscPrintf(PETSCFEM_COMM_WORLD,
                            "fails to read integer per-element"
                            " property %d,\n at line \"%s\"",
                            jprop+1,line);
                PFEMERRQ("\n");
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
        DONE:;
        } else {
          AutoString datas;
          datas.set(data).deblank();
          data = datas.str();
          TGETOPTDEF(thash,int,use_hdf5,0);
          Autobuf *tempo = abuf_create();
          FileStack *file_connect=NULL;
          if (!myrank) {
            file_connect = new FileStack(data);
            ierro = !file_connect->ok();
          }
          CHECK_PAR_ERR(ierro,"Error opening connectivity file");
          if (!myrank) {
            nelem=0;
            while (!file_connect->get_line(line)) {
              // reading element connectivities
              for (int jel=0; jel<nel; jel++) {
                token =  strtok(( jel==0 ? line : NULL),bsp);
                if (!token) file_connect->print();
                PETSCFEM_ASSERT(token,
                                "Error reading element connectivities at\n"
                                "%s:%d: \"%s\"",file_connect->file_name(),
                                file_connect->line_number(),
                                file_connect->line_read());
                sscanf(token ,"%d",&node);
                icorow[jel]= node;
              }

              // reading element properties
              for (int jprop=0; jprop<nelprops; jprop++) {
                token = strtok(NULL,bsp);
                if (token==NULL) {
                  PetscPrintf(PETSCFEM_COMM_WORLD,
                              "fails to read per-element property %d,\n at line \"%s\"",
                              jprop+1,line);
                  PFEMERRQ("\n");
                }
                sscanf(token,"%lf",proprow+jprop);
              }

              // reading integer element properties
              for (int jprop=0; jprop<neliprops; jprop++) {
                token = strtok(NULL,bsp);
                if (token==NULL) {
                  PetscPrintf(PETSCFEM_COMM_WORLD,
                              "fails to read integer per-element"
                              " property %d,\n at line \"%s\"",
                              jprop+1,line);
                  PFEMERRQ("\n");
                }
                sscanf(token,"%d",iproprow+jprop);
              }
	
              // Copying to buffer
              ierr = abuf_cat_buf (tempo,(unsigned char *)icorow,nel*sizeof(int));
              CHKERRQ(ierr);
              ierr = abuf_cat_buf (tempo,(unsigned char *)proprow,nelprops*sizeof(double));
              CHKERRQ(ierr);
              ierr = abuf_cat_buf (tempo,(unsigned char *)iproprow,neliprops*sizeof(int));
              CHKERRQ(ierr);
              nelem++;
            }
            ierro = file_connect->last_error()!=FileStack::eof;
            if (ierro) printf("Couldn't process correctly connectivity file%s\n",
                              file_connect->file_name());
            file_connect->close();
            delete file_connect;
          }
          CHECK_PAR_ERR(ierro,"Error reading connectivity file");
          ierr = MPI_Bcast (&nelem,1,MPI_INT,0,PETSCFEM_COMM_WORLD);
          // Resize `da_icone' in other processors
          if (myrank) ierr =  abuf_set (tempo,nelem*rowsize,0);
	  
          // Broadcast all `icone' data using MPI
          ierr = MPI_Bcast (abuf_data(tempo),nelem*rowsize,MPI_CHAR,0,PETSCFEM_COMM_WORLD);

          // This should be done AFTER reading the nodes 
          // Set all nodes that are connected to an element as degrees of freedom
          for (int e=0; e<nelem; e++) {
            memcpy (icorow,abuf_data(tempo)+e*rowsize,nel*sizeof(int));
            for (int jel=0; jel<nel; jel++) {
              node = icorow[jel];
              for (int kdof=1; kdof<=ndof; kdof++) {
                edof = dofmap->edof(node,kdof);
                dofmap->id->set_elem(edof,edof,1.);
              }
            }
            // Copying data from tempo tu `da_icone'
            int indxi = da_append(da_icone,abuf_data(tempo)+e*rowsize);
            if ( indxi==-1 ) PFEMERRQ("Insufficient memory reading elements");
          }
          abuf_destroy(tempo);
        }	
    
        elemsetnum++;
        delete[] icorow;

        nelem = da_length(da_icone);

        elemprops = new double[nelem*nelprops];
        elemiprops = new int[nelem*neliprops];
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
      }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
      // Bless with the appropriate type
      // adding a line CHECK_ELEMSET_TYPE(type) for each
      // type of elemset to be added
      bless_elemset(type,elemset);
      if(!elemset) PETSCFEM_ERROR("not known elemset type: \"%s\"\n",type);

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

      string name = Elemset::anon;
      ::get_string(thash,"name",name,1);
      elemset->register_name(name,type);
      thash->register_name(elemset->name());
      if (myrank==0) thash->print("Table of properties:");

      elemset->read(fstack);
      elemset->initialize();

      PetscPrintf(PETSCFEM_COMM_WORLD,
		  "elemset number %d, name \"%s\", pointer %p, number of elements %d\n",
		  elemsetnum,elemset->name(),(void*)elemset,nelem);
      
      // Append to the list
      da_append(mesh->elemsetlist,&elemset);
      PetscPrintf(PETSCFEM_COMM_WORLD,"Ends reading  elemset\n");

      
    } else if (!strcmp(token,"end_elemsets")) {

      // nothing is done here
      PetscPrintf(PETSCFEM_COMM_WORLD,"End elemsets section\n");

    } else if (!strcmp(token,"fixa")) {

      PetscPrintf(PETSCFEM_COMM_WORLD," -- Reading fixations\n"); 
      // Read fixations
      // dofmap->fixa = da_create(sizeof(fixa_entry));
      // fixa_entry fe;
      nfixa=0;
      while (1) {
	ierr = fstack->get_line(line);
	PETSCFEM_ASSERT0(ierr==0,"Can't find __END_FIXA__ tag");  
	astr_copy_s(linecopy, line);
	if (strstr(line,"__END_FIXA__")) break;
	nfixa++;
	int nr=sscanf(line,"%d %d %lf",&node,&kdof,&dval);
	if (nr !=3) {
	  PetscPrintf(PETSCFEM_COMM_WORLD,
		      "Error reading fixations, for fixation %d\nline: %s\n",
		      nfixa,astr_chars(linecopy));
	  CHKERRQ(1);
	}

	dofmap->get_row(node,kdof,row);
	if (row.size()!=1) {
	  PetscPrintf(PETSCFEM_COMM_WORLD,
		      "In line: %s\nFixation %d, imposed on an invalid"
		      " node/field combination.\n",
		      astr_chars(linecopy),nfixa);
	}
	  
	int keq = row.begin()->first;
	PETSCFEM_ASSERT(row.begin()->second == 1.,
			"Fixation imposed on a bad node/field combination\n"
			"node: %d, field: %d\n"
			"%s:%d: \"%s\"",
			node,kdof,
			fstack->file_name(),
			fstack->line_number(),
			fstack->line_read());
	
	map<int,int>::iterator q = dofmap->fixed_dofs.find(keq);
	if (q!=dofmap->fixed_dofs.end()) {
	  int j = q->second;
	  // delete dofmap->fixed[j];
	  dofmap->fixed[j] = fixation_entry(dval);
	} else {
	  dofmap->fixed.push_back(fixation_entry(dval));
	  dofmap->fixed_dofs[keq]=dofmap->fixed.size()-1;
	}
      }
      PetscPrintf(PETSCFEM_COMM_WORLD,"Total fixations: %d\n",nfixa);

    } else if (!strcmp(token,"fixa_amplitude")) {

      PetscPrintf(PETSCFEM_COMM_WORLD," -- Reading fixa_amplitude section\n"); 
      // next token is the identifier of the amplitude function
      token = strtok(NULL,bsp);
      PETSCFEM_ASSERT0(token,"Couldn't find name for fixa_amplitude section");  
      char *label = new char[strlen(token)+1];
      strcpy(label,token);

      Amplitude *amp = Amplitude::old_factory(label,*fstack);

      string datah5,dvalsh5;
      int use_hdf5_fixa=0;
      if (!amp) {
	// Read options table
	TextHashTable *thash = new TextHashTable;
	read_hash_table(fstack,thash);
        TGETOPTDEF_S(thash,string,data,NONE);
        datah5 = data;
        TGETOPTDEF(thash,int,use_hdf5,-1);
        use_hdf5_fixa = use_hdf5;
        TGETOPTDEF_S(thash,string,datavals,NONE);
        dvalsh5 = datavals;
            
	// create a new Amplitude
	amp = Amplitude::factory(label,thash);
      }
      // amplitude_list:= stores a list of all the amplitudes defined 
      amplitude_list.push_back(amp);

      TGETOPTDEF(thash,int,use_hdf5,0);
      TGETOPTDEF(thash,int,use_hdf5_auto,0);
      int hdf5 = use_hdf5;
      if (use_hdf5_auto && ishdf5(datah5.c_str())) hdf5=1;
      if (use_hdf5_fixa>=0) hdf5 = use_hdf5_fixa;
      if (hdf5) {
        const char *data = datah5.c_str();
        dvector<int> dvfixa;
        dvector<double> dvvals;
        if (!myrank) {
          h5_dvector_read(data,dvfixa);
          dvfixa.defrag();
          nfixa = dvfixa.size();
          PETSCFEM_ASSERT(nfixa%2==0,
                          "bad elemset HDF5 dataset size, "
                          "size %d",nfixa);
          nfixa /= 2;
          dvfixa.reshape(2,nfixa,2);
          printf("data %s, nfixa %d\n",data,nfixa);

          PETSCFEM_ASSERT0(dvalsh5!="NONE",
                           "datavals entry is required if use_hdf5");  
          h5_dvector_read(dvalsh5.c_str(),dvvals);
          PETSCFEM_ASSERT(dvvals.size()==nfixa,
                          "datavals must have the same elemens "
                          "as fixa entries. "
                          "fixa entries %d, datavals entries %d",
                          nfixa,dvvals.size());  
        }
        ierr = MPI_Bcast(&nfixa,1,MPI_INT,0,PETSCFEM_COMM_WORLD);
        dvector_clone_parallel(dvfixa);
        dvfixa.defrag();

        dvector_clone_parallel(dvvals);
        dvvals.defrag();
        PetscPrintf(PETSCFEM_COMM_WORLD,
                    "Reading fixations HDF5. %d fixations read\n",nfixa);
        exit(0);

#if 0
        dofmap->get_row(node,kdof,row);
        if (row.size()!=1) {
          PetscPrintf(PETSCFEM_COMM_WORLD,
                      "In line: %s\nFixation %d, imposed on an"
                      " invalid node/field combination.\n",
                      astr_chars(linecopy),nfixa);
        }
	  
        int keq = row.begin()->first;
        PETSCFEM_ASSERT(row.begin()->second == 1.,
                        "Fixation imposed on a bad node/field combination\n"
                        "node: %d, field: %d\n"
                        "%s:%d: \"%s\"",
                        node,kdof,
                        fstack->file_name(),
                        fstack->line_number(),
                        fstack->line_read());
	
        edof = dofmap->edof(node,kdof);
        map<int,int>::iterator q = dofmap->fixed_dofs.find(keq);
        if (q!=dofmap->fixed_dofs.end()) {
          int j = q->second;
          // delete dofmap->fixed[j];
          dofmap->fixed[j] = fixation_entry(dval,amp,edof);
        } else {
          dofmap->fixed.push_back(fixation_entry(dval,amp,edof));
          dofmap->fixed_dofs[keq]=dofmap->fixed.size()-1;
        }
#endif        
      } else {
        nfixa=0;
        while (1) {
          ierr = fstack->get_line(line);
          PETSCFEM_ASSERT0(ierr==0,"Can't find __END_FIXA__ tag");  
          astr_copy_s(linecopy, line);
          if (strstr(line,"__END_FIXA__")) break;
          nfixa++;
          int nr=sscanf(line,"%d %d %lf",&node,&kdof,&dval);
          if (nr !=3) {
            PetscPrintf(PETSCFEM_COMM_WORLD,
                        "Error reading fixations, for fixation %d\n",nfixa);
            CHKERRQ(1);
          }

          dofmap->get_row(node,kdof,row);
          if (row.size()!=1) {
            PetscPrintf(PETSCFEM_COMM_WORLD,
                        "In line: %s\nFixation %d, imposed on an"
                        " invalid node/field combination.\n",
                        astr_chars(linecopy),nfixa);
          }
	  
          int keq = row.begin()->first;
          PETSCFEM_ASSERT(row.begin()->second == 1.,
                          "Fixation imposed on a bad node/field combination\n"
                          "node: %d, field: %d\n"
                          "%s:%d: \"%s\"",
                          node,kdof,
                          fstack->file_name(),
                          fstack->line_number(),
                          fstack->line_read());
	
          edof = dofmap->edof(node,kdof);
          map<int,int>::iterator q = dofmap->fixed_dofs.find(keq);
          if (q!=dofmap->fixed_dofs.end()) {
            int j = q->second;
            // delete dofmap->fixed[j];
            dofmap->fixed[j] = fixation_entry(dval,amp,edof);
          } else {
            dofmap->fixed.push_back(fixation_entry(dval,amp,edof));
            dofmap->fixed_dofs[keq]=dofmap->fixed.size()-1;
          }
        }
        PetscPrintf(PETSCFEM_COMM_WORLD,
                    "Total fixations with temporal amplitude: %d\n",nfixa);
      }

    } else if (!strcmp(token,"constraint")) {

      PetscPrintf(PETSCFEM_COMM_WORLD," -- Reading constraint section\n"); 
#define ERRLINE								\
      if (ierr==-1) {							\
         PetscPrintf(PETSCFEM_COMM_WORLD,"Error reading line:\n\"%s\"\n",	\
		     astr_chars(linecopy));				\
         CHKERRQ(1);							\
      }

      // Read constraints
      int node1,kdof1,nconstr=0,nlindep=0;
      Constraint constraint;
      while (1) {
	
	ierr = fstack->get_line(line);
	astr_copy_s(linecopy, line);
	PETSCFEM_ASSERT0(ierr==0,"Can't find __END_CONSTRAINT__ tag");  
	if (strstr(line,"__END_CONSTRAINT__")) break;
	nconstr++;
	constraint.empty();
	double coef;
	int node,field;
	rflag=0;
// #define PETSCFEM_WARNING(bool_cond,templ,...) 
// if (!(bool_cond)) { PetscPrintf(PETSCFEM_COMM_WORLD, 
// 				"PETSc-FEM warning: " templ,__VA_ARGS__);}
	while (1) { 
	  ierr = readval(rflag,line,coef); ERRLINE;
	  if (!ierr) break;
	  ierr = readval(rflag,line,node); ERRLINE;
	  ierr = readval(rflag,line,field); ERRLINE;
	  PETSCFEM_ASSERT(node<=nnod,"read_mesh: "
			  "Read node= %d greater than nnod= %d\n"
			  "%s:%d: \"%s\"",node,nnod,fstack->file_name(),
			  fstack->line_number(),
			  fstack->line_read());
	  PETSCFEM_ASSERT(field<=ndof,"read_mesh: "
			  "Read field= %d greater than ndof= %d\n"
			  "%s:%d: \"%s\"",field,ndof,fstack->file_name(),
			  fstack->line_number(),
			  fstack->line_read());
	  constraint.add_entry(node,field,coef);
	}
	int return_status = dofmap->set_constraint(constraint);
	if (return_status==1) {
	  PetscPrintf(PETSCFEM_COMM_WORLD,
		      "Linearly dependent fixation discarded. \n"
		      "%s:%d: \"%s\"\n",
		      fstack->file_name(),
		      fstack->line_number(),
		      fstack->line_read());
	  nlindep++;
	} else if (return_status==2) {
	  PETSCFEM_ERROR("Couldn't find a free dof in constraint list.\n"
			 "%s:%d: \"%s\"\n",
			 fstack->file_name(),
			 fstack->line_number(),
			 fstack->line_read());
	}
      }
      // dofmap->id->print("");
      PetscPrintf(PETSCFEM_COMM_WORLD,"Total entered constraint lines: %d.\n",nconstr);
      PetscPrintf(PETSCFEM_COMM_WORLD,"-- Linearly dependent %d, linearly independent %d\n",
		  nlindep,nconstr-nlindep);
#undef ERRLINE

    } else {

      // These are objects that read themselves (OOP)
      // In a future all objects will be read this way and
      // will provide, unknowns, fields, equations and constraints. 
      string basic_object_type(token);
      BasicObject *obj = BasicObject::factory(basic_object_type);
      if (obj) {
	obj->read(fstack,mesh,dofmap);
      } else {
        PETSCFEM_ERROR("Bad section name in data file.\n"
		    "line: \"%s\"\n",line);  
      }
    }

  }

  PETSCFEM_ASSERT(fstack->last_error()==FileStack::eof,
		  "Couldn't process correctly main data "
		  "file \"%s\"\n",fstack->file_name());  
  PETSCFEM_ASSERT0(nnod>0,"No nodes read");  
  PETSCFEM_ASSERT0(da_length(mesh->elemsetlist)>0,
		   "No elemsets read");  
    
  fstack->close();
  delete fstack;

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
  PetscPrintf(PETSCFEM_COMM_WORLD,"size: %d\n",size);
  if (size==1 || proc_weights == NULL) {

    // If there is not  a processor weight table, then 
    for (int proc=0; proc<size; proc++) {
      tpwgts[proc] = 1./float(size);
    }

  } else {

    weights_file = new FileStack(proc_weights);
    ierro = !weights_file->ok();
    if (!ierro) {
      weights_file->quiet = !myrank;
      float sumw=0.;
      for (int proc=0; proc<size; proc++) {
	ierr = weights_file->get_line(line);
	PETSCFEM_ASSERT(ierr==0,"Can't find line for proc %d in file %s",
			proc,proc_weights);  
	sscanf(line,"%f",&tpwgts[proc]);
	PETSCFEM_ASSERT0(tpwgts[proc]>=0.,
			 "Processor weight must be >= 0.");  
	sumw += tpwgts[proc];
      }
      PetscPrintf(PETSCFEM_COMM_WORLD,"total weight: %f\n",sumw);
      PETSCFEM_ASSERT0(sumw>0.,"Total processor weight must be > 0.");  
      for (int proc=0; proc<size; proc++) {
	tpwgts[proc] /= sumw;
	PetscPrintf(PETSCFEM_COMM_WORLD," proc: %d, w: %f\n",proc,tpwgts[proc]);
      }
      weights_file->close();
      delete weights_file;
    }
    CHECK_PAR_ERR(ierro,"Error reading weights file. ");
    
  }
    
  //o Print hostnames for nodes participating in this run
  TGETOPTDEF(mesh->global_options,int,print_hostnames,0);
  if (print_hostnames) {
#define MAX_NAME_LEN 200
    char line[MAX_NAME_LEN];
    ierr = gethostname(line,MAX_NAME_LEN);
    if (ierr) sprintf(line,"<unknown>");
    KeyedOutputBuffer buff;
    buff.cat_printf("[%d] hostname %s",myrank,line);
    buff.push(myrank);
    buff.flush();
  }

  // nelemsets:= total number of elemsets in the mesh
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
  PFEMERRCQ(size>1 && !is_any_fat,
	    "At least one elemset has to be set as fat (is_fat property)\n"
	    "if number of processors is np>1\n");

  // cumulated sum 
  for (node=0; node<nnod; node++) 
    n2eptr[node+1] += n2eptr[node];
  // PetscPrintf(PETSCFEM_COMM_WORLD,"n2eptr[nnod]: %d\n",n2eptr[nnod]);

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

#if 0
  for (node=0; node<nnod; node++) {
    printf("node %d: ",node+1);
    for (int jpos=n2eptr[node]; jpos<n2eptr[node+1]; jpos++) 
      printf(" %d",node2elem[jpos]);
    printf("\n");
  }
#endif
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
  
//    // o Do not coalesce elements if the number of elements if below this
//    // limit. 
//    TGE TOPTDEF(mesh->global_options,int,min_partgraph_vertices,0);
//    // o Coalesce elements by this ratio if number of elements 
//    TGE TOPTDEF(mesh->global_options,int,partgraph_vert_to_nodes_ratio,INF);

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
    PETSCFEM_ERROR("partitioning method not known \"%s\"\n",
		   partitioning_method.c_str());
  }
  PetscPrintf(PETSCFEM_COMM_WORLD,"Starts partitioning.\n"); 

    
  int *vpart = new int[nelemfat];
  if (partflag==0 || partflag==2) {

    metis_part(nelemfat,mesh,nelemsets,vpart,
	       nelemsetptr,n2eptr,node2elem,size,myrank,
	       partflag,tpwgts,max_partgraph_vertices,
	       iisd_subpart,print_partitioning_statistics);

  } else if (partflag==1) {

    int *mnel = new int[size+1];
    if (size>1) {
      if (myrank==0) printf("Hitchicking partition, partflag = %d\n",partflag);
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
    if (myrank==0) {
      for (int j=0; j < nelemfat; j++) {
	vpart[j] = rand() % size;
	if (partflag==4)
	  vpart[j] = int(double(j)/double(nelemfat)*double(size));
	// Just in case random functions are too close to the limits
	// (In theory, rand() should not touch the limits 0, RAND_MAX)
	if (vpart[j]>=size) vpart[j]=size-1;
	if (vpart[j]<0) vpart[j]=0;
	// printf("element %d in [%d]\n",j+1,vpart[j]);
      }
    }
    ierr = MPI_Bcast (vpart,nelemfat,MPI_INT,0,PETSCFEM_COMM_WORLD);
    
  } else assert(0); // something went wrong

  PetscPrintf(PETSCFEM_COMM_WORLD,"Ends partitioning.\n");

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

  dvector<double> ii_stat, ii_stat2;
  ii_stat.a_resize(2,size,size);
  ii_stat2.a_resize(1,size);
  if (print_nodal_partitioning)
    PetscPrintf(PETSCFEM_COMM_WORLD,
                "\nNodal partitioning (node/processor): \n");

  // Node interface between processor statistics
  int c1,c2,D1,D2,P1,P2;
  ii_stat.set(0.0);
  for (node=0; node<nnod; node++) {
    ii_stat2.set(0.0);
    for (c1 = n2eptr[node]; c1 < n2eptr[node+1]; c1++) {
      D1 = vpart[node2elem[c1]];
      ii_stat2.e(D1) += 1.0;
    }
    double mass = n2eptr[node+1]-n2eptr[node];
    if (mass!=0.0) 
      for (int p=0; p<size; p++) ii_stat2.e(p) /= mass;
    for (D1=0; D1<size; D1++)
      for (D2=0; D2<size; D2++) 
	ii_stat.e(D1,D2) += ii_stat2.e(D1)*ii_stat2.e(D2);
  }

  //o For multi-core architectures this represent the number
  //  of cores per node. 
  GETOPTDEF(int,ncore,0);
  PETSCFEM_ASSERT(ncore>=0,"ncore must be non-negative, given %d",
                  ncore);

  GETOPTDEF(int,multicore_rand_part,0);

  dvector<int> procmap;
  procmap.mono(size).reshape(1,size);
  for (int j=0; j<size; j++) procmap.e(j) = j;

  if (size>1 && ncore>0) {

    PETSCFEM_ASSERT(size%ncore==0,
                    "if ncore is given, numer of "
                    "processor `size' must be a multiple of ncore.\n"
                    "size %d, ncore %d",size,ncore);

    if (!myrank) {

      // Using graph-matching algorithm
      dvector<double> bw;
      // Bandwidths for connection inside and outside the nodes
      double bw1 = 100.0, bw2 = 1.0;
      bw.mono(size*size).reshape(2,size,size).set(0.0);
      for (int j=0; j<size; j++) {
        int nodej = j/ncore;
        for (int k=0; k<size; k++) {
        int nodek = k/ncore;
        if (j == k) continue;
        bw.e(j,k) = (nodej==nodek? bw1 : bw2);
        }
      }

      // Compute statistics without matching algorithm
      double pmax=NAN, psum=NAN;
      perfo(bw,ii_stat,procmap,pmax,psum);
      printf("natural partition: max %g, sum %g\n",
             pmax,psum);

#if 1
      // Apply matching graph algorithm
      match_graph(bw,ii_stat,procmap);
      perfo(bw,ii_stat,procmap,pmax,psum);
      printf("match graph algorithm: max %g, sum %g\n",
             pmax,psum);
#endif

#if 1
      // Random, just to verify
      dvector<int> procmap_rand;
      procmap_rand.clone(procmap);
      procmap_rand.defrag();
      int *proc_p = procmap_rand.buff();
      random_shuffle(proc_p,proc_p+size);
      perfo(bw,ii_stat,procmap_rand,pmax,psum);
      printf("random permut: max %g, sum %g\n",pmax,psum);
#endif

#if 1
      // Apply partitioning specific for tree like graph
      repart(ii_stat,procmap,ncore);
      perfo(bw,ii_stat,procmap,pmax,psum);
      printf("tree partitioning: max %g, sum %g\n",
             pmax,psum);
#endif

      if (multicore_rand_part) {
        printf("Using random allocation for multicore partitioning...\n");
        procmap.clone(procmap_rand);
      }

#if 0
      for (int j=0; j<size; j++)
        printf("domain %d, sent to processor %d\n",
               j,procmap.e(j));
#endif

#if 0
      printf("Bandwidth load performance [w/o algorithm, with, "
             "improvement ratio]\n");
      printf("  bw ratio MAX: [%g, %g, %g]\n",
             pmax_mg0,pmax_mg1,pmax_mg1/pmax_mg0);
      printf("  bw ratio SUM: [%g, %g, %g]\n",
             psum_mg0,psum_mg1,psum_mg1/psum_mg0);
#endif
      for (int j=0; j<nelemfat; j++) {
        assert(vpart[j]>=0 && vpart[j]<size);
        vpart[j] = procmap.e(vpart[j]);
      }
    }
 
    ierr = MPI_Bcast(vpart,nelemfat,MPI_INT,0,PETSCFEM_COMM_WORLD);
    if (!myrank) {
      dvector<double> ii_stat_cpy;
      ii_stat_cpy.clone(ii_stat);
      for (int dj=0; dj<size; dj++) {
        int pj = procmap.e(dj);
        for (int dk=0; dk<size; dk++) {
          int pk = procmap.e(dk);
          ii_stat_cpy.e(pj,pk) = ii_stat.e(dj,dk);
        }
      }
      ii_stat.clone(ii_stat_cpy);
    }
    ierr = MPI_Bcast (ii_stat.buff(),size*size,
                      MPI_INT,0,PETSCFEM_COMM_WORLD);
  } 

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

  if (myrank == 0 && size > 1) {
    PetscPrintf(PETSCFEM_COMM_WORLD,
                "---\nInter-processor node connections\n");
    // compute inverse map
    dvector<int> dommap;
    dommap.clone(procmap);
    for (int dj=0; dj<size; dj++) {
      int pj = procmap.e(dj);
      dommap.e(pj) = dj;
    }
    for (P1=0; P1<size; P1++) {
      D1 = dommap.e(P1);
      for (P2=0; P2<size; P2++) {
        D2 = dommap.e(P2);
	printf("[%d]-[%d] %.2f\n",D1,D2,ii_stat.e(D1,D2));
      }
    }
    PetscPrintf(PETSCFEM_COMM_WORLD,"\n");
  }
  ii_stat.clear();
  ii_stat2.clear();

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
      PetscPrintf(PETSCFEM_COMM_WORLD,
		  "%d   %d\n",node+1,npart[node]);
  }
  if (print_nodal_partitioning)
    PetscPrintf(PETSCFEM_COMM_WORLD,"End nodal partitioning table\n\n");

  if (node_not_connected_to_fat)
    PetscPrintf(PETSCFEM_COMM_WORLD,"warning! there are %d "
		"nodes not linked to any \"fat\" elemset. \n"
		"This induces artificial numbering. [But may be OK]\n",
		node_not_connected_to_fat);

  //o Prints element partitioning. 
  GETOPTDEF(int,debug_element_partitioning,0);

  // Define the eparts of each fat elemset as the corresponding part
  // of vpart. 
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
      PetscPrintf(PETSCFEM_COMM_WORLD,
		  "Elemset \"%s\", type \"%s\", nelem %d\n."
		  "Partitioning table: \n",elemset->name(),elemset->type,elemset->nelem);
      for (int kk=0; kk<elemset->nelem; kk++) 
	PetscPrintf(PETSCFEM_COMM_WORLD,"%d -> %d\n",kk,epart[kk]);
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

  // Number degrees of freedom
  int proc,*neqproc, *startproc;
  neqproc= new int[size+1];
  startproc= new int[size+1];
  int *perm = new int [ndof*nnod];
  for (k=0; k<ndof*nnod; k++) {
    perm[k]=0;
  }


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
  PetscPrintf(PETSCFEM_COMM_WORLD,
	      "Total number of degrees of freedom neq:     %d\n"
	      "Total number of independent fixations neqf: %d\n",
	      neq,dofmap->neqf);

  if (size>1) {
    for (proc=0; proc<size; proc++) 
      PetscPrintf(PETSCFEM_COMM_WORLD,
		  "[%d] number of dof's: %d\n",proc,neqproc[proc]);
  }
  

#if 0
  // This is the old version related to the remap_cols() bug. A leak
  // of memory due to the memory management of the STL map routines. 
  dofmap->id->remap_cols(perm);
#else
  idmap *idnew = new idmap(nnod*ndof,NULL_MAP);
  dofmap->id->remap_cols(perm,*idnew);
  DELETE_SCLR(dofmap->id);
  dofmap->id = idnew;
#endif
 

  //o Checks that the  #idmap#  has been correctly generated. 
  TGETOPTDEF(mesh->global_options,int,check_dofmap_id,0);
  if (check_dofmap_id && myrank==0) {
    dofmap->id->check();
    dofmap->id->print_statistics();
  }


  //o Prints the dofmap  #idmap#  object. 
  TGETOPTDEF(mesh->global_options,int,print_dofmap_id,0);
  if (print_dofmap_id && myrank==0) {
    dofmap->id->print("dofmap->id: \n");
  }

  dofmap->freeze();

  map<int,int>::iterator jj;

#if 0  // debug:=
  {
    map<int,int>::iterator jjjj;
    for (jjjj=dofmap->fixed_dofs.begin(); jjjj!=end_fix; jjjj++) 
      printf("%d -> fix %d\n",jjjj->first,jjjj->second);

    for (int kkk=0; kkk<nnod*ndof; kkk++) {
      printf("%d -> perm[%d]\n",kkk,perm[kkk]);
    }
  
  }
#endif

  int nfixed = dofmap->fixed.size();
  vector<fixation_entry> fixed_remapped(nfixed);
  //  fixed_remapped.reserve(nfixed);
  for (jj=dofmap->fixed_dofs.begin(); jj!=dofmap->fixed_dofs.end(); jj++) {
    int newjeq = perm[jj->first -1];
    if (newjeq==0) continue;
    int indx= newjeq - neq - 1;
    assert(0<=indx && indx<nfixed); //This should be in this range. Otherwise it falls
				    // off the `fixed_remapped' vector.
    fixed_remapped[indx] =
      (dofmap->fixed)[jj->second];
  }


  // swap fixed with fixed_remapped (reordered)
  dofmap->fixed.swap(fixed_remapped);
  VOID_IT(fixed_remapped);
  VOID_IT(dofmap->fixed_dofs);

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
      elemset->epart = new int[nelem];
      const int *conn; 
      for (iele=0; iele<nelem; iele++) {
	// Decide in which processor will be computed this element
        int nell = elemset->real_nodes(iele,conn);
        PETSCFEM_ASSERT0(nell>0,"Elements should have at least one node");  
	node = conn[0];
	proc = npart[node-1];
	if (proc<1 || proc>size) {
	  PetscPrintf(PETSCFEM_COMM_WORLD,
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
    { // Make a local scope.
      // Compute the epart2 and epart_p vectors. 
      // Later epart2 could be avoided and only using epart.
      // Right now we leave both, but it is coded so that
      // later we could simply assign epart2=epart and it should work. 
      vector<int> &epart_p = elemset->epart_p;
      epart_p.resize(size+1);
      int *epart2 = new int[elemset->nelem];
      elemset->epart2 = epart2;
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
      elemset->e1 = epart_p[myrank];
      elemset->e2 = epart_p[myrank+1];
#if 0
      if (!myrank) {
	for (int k=0; k<elemset->nelem; k++) {
	  printf("elem %d, epart %d, epart2 %d\n",k,epart[k],epart2[k]);
	}
      }
      PetscFinalize();
      exit(0);
#endif    
    }
    // ghost_elems:= These are elements that have related dof's on the
    // processor, but they don't live on the processor.
    // Loops for defining profiles of matrices must loop over them
    // also. 

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
      elemset->local_store = new void *[elemset->nelem_here];
      for (int j=0; j<elemset->nelem_here; j++) {
	elemset->local_store[j]=NULL;
      }
    }

    
    // da_sort (elemset->ghost_elems,int_cmp,NULL);
    int nghostel = da_length(elemset->ghost_elems);
    if (nghostel>0) {
      int *dap = (int *)da_ref(elemset->ghost_elems,0);
      sort(dap,dap+nghostel);
    }

    if (size>1) {
      PetscPrintf(PETSCFEM_COMM_WORLD,
		  "For elemset type \"%s\", name \"%s\"\n",
		  elemset->type,elemset->name());
      PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,
			      "On processor [%d], %d local elements,"
			      " %d ghost elements.\n",
			      myrank,elemset->nelem_here,nghostel);
      PetscSynchronizedFlush(PETSCFEM_COMM_WORLD);
    }

  }
  
  // Convert ghost_dof_set en dofmap
  set<int>::iterator it;
  for (it=ghost_dof_set.begin(); it!=ghost_dof_set.end(); it++) 
    dofmap->ghost_dofs.push_back(*it);
  PetscSynchronizedFlush(PETSCFEM_COMM_WORLD);
  VOID_IT(ghost_dof_set);
  int nghost_dofs = dofmap->ghost_dofs.size();

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

  // Defines certain quantities in dofmap
  dofmap->neq= neq;
  dofmap->startproc = startproc;
  dofmap->neqproc = neqproc;
  dofmap->size = size;
  dofmap->npart = npart;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  //                        DEFINE SCATTER
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

  PetscPrintf(PETSCFEM_COMM_WORLD,"Defining scatters...\n");

  // Create MPI vectors
  //  Vec x,xseq;
  Vec x,ghost_vec,xseq;
  IS is_ghost_glob,is_ghost_loc,is_print;

  dofmap->create_MPI_vector(x);
  dofmap->create_MPI_ghost_vector(ghost_vec);

  int nghost_tot;
  ierr = VecGetSize(ghost_vec,&nghost_tot); CHKERRQ(ierr);

  ierr = ISCreateGeneral(PETSCFEM_COMM_WORLD,nghost_dofs,
			 dofmap->ghost_dofs.data(),
                         PETSC_COPY_VALUES,
			 &is_ghost_glob);  CHKERRQ(ierr); 
  ierr = ISCreateStride(PETSCFEM_COMM_WORLD,nghost_dofs,0,1,&is_ghost_loc);

  ierr = VecScatterCreate(x,is_ghost_glob,ghost_vec,is_ghost_loc,
			  &dofmap->ghost_scatter); CHKERRQ(ierr); 

  ierr = ISDestroy(&is_ghost_glob); CHKERRQ(ierr); 
  ierr = ISDestroy(&is_ghost_loc); CHKERRQ(ierr); 

  int neql = (myrank==0 ? dofmap->neq : 0);
  ierr = VecCreateSeq(PETSC_COMM_SELF,neql,&xseq);  CHKERRQ(ierr);
  
  ierr = ISCreateStride(PETSCFEM_COMM_WORLD,neql,0,1,&is_print);
  ierr = VecScatterCreate(x,is_print,xseq,is_print,
			  &dofmap->scatter_print); CHKERRQ(ierr); 
  
  ierr = ISDestroy(&is_print); CHKERRQ(ierr); 

#if 0
  for (int jj=dof1; jj<=dof2; jj++) {
    VecSetValue(x,jj-1,double(jj),INSERT_VALUES);
  }
  PetscPrintf(PETSCFEM_COMM_WORLD,"vector x en read_mesh\n");
  ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD); CHKERRA(ierr);
  PetscPrintf(PETSCFEM_COMM_WORLD,"================\n");
  ierr = VecScatterBegin(dofmap->ghost_scatter,x,ghost_vec,
			 INSERT_VALUES,SCATTER_FORWARD); CHKERRA(ierr); 
  ierr = VecScatterEnd(dofmap->ghost_scatter,x,ghost_vec,
		       INSERT_VALUES,SCATTER_FORWARD); CHKERRA(ierr); 
  
  double *array;
  ierr = VecGetArray(ghost_vec,&array); CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,"Local ghost values on [%d]\n",myrank);
  for (int jj=0; jj<nghost_dofs; jj++ ) 
    PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,"local %d, dof %d  -> %g\n",
			    jj,dofmap->ghost_dofs[jj],array[jj]);
  PetscSynchronizedFlush(PETSCFEM_COMM_WORLD);
#endif

  ierr = VecDestroy(&x);
  ierr = VecDestroy(&ghost_vec);
  ierr = VecDestroy(&xseq);

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

  PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,"On processor [%d], local values\n",myrank);
  for (int jj=0; jj<neqproc[myrank]; jj++ ) {
    PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,"local %d, global %d -> %g\n",
			    jj,startproc[myrank]+jj,array[jj]);
  }

  PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,"Local ghost values\n");
  for (int jj=0; jj<nghost_dofs; jj++ ) {
    int local = neqproc[myrank]+jj;
    PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,"local %d, global %d -> %g\n",
			    local,dofmap->ghost_dofs[jj],array[local]);
  }

  PetscSynchronizedFlush(PETSCFEM_COMM_WORLD);
  ierr = VecRestoreArray(lx,&array);CHKERRQ(ierr);
  ierr = VecGhostRestoreLocalForm(gx,&lx);CHKERRQ(ierr); 
  PetscFinalize();
  exit(0);
#endif 

#if 0
  // para debuggear
  ierr = VecScatterBegin(scatter,x,xseq,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr); 
  ierr = VecScatterEnd(scatter,x,xseq,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr); 
  PetscPrintf(PETSCFEM_COMM_WORLD,"Despues del scatter\n");


  // debug
  ierr = VecGetArray(xseq,&sstate); CHKERRQ(ierr); 
  PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,"On processor %d ------------------\n",
			  myrank+1);
  for (int k=0; k<neq; k++) {
    PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,"%d -> %f\n",k+1,sstate[k]);
  }
  PetscSynchronizedFlush(PETSCFEM_COMM_WORLD);
  ierr = VecRestoreArray(xseq,&sstate); CHKERRQ(ierr); 

  PetscFinalize();
  exit(0);
#endif  

#if 0
  PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,
			  "On processor %d \n",myrank+1);
  for (k=0; k<ndofhere; k++) {
    PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,
			    "%d -> %d\n",k+1,dof_here_list[k]);
  }
  PetscSynchronizedFlush(PETSCFEM_COMM_WORLD);
#endif
  
  PetscPrintf(PETSCFEM_COMM_WORLD,"Total number of dof's: %d\n",neq);

#if 0
  int nloads=0;
  fstack->get_line(line);
  if (strncmp("loads",line,5) ) {
    PetscPrintf(PETSCFEM_COMM_WORLD,"Couldn't find <fixa> tag line\n");
    exit(1);
  }
  while (1) {
    fstack->get_line(line);
    if (strstr(line,"__END_LOADS__")) break;
    nloads++;
    int nr=sscanf(line,"%d %d %lf",&node,&kdof,&dval);
    if (nr !=3) {
      PetscPrintf(PETSCFEM_COMM_WORLD,"Error reading LOADS, for loads %d\n",nfixa);
    }
    //LOAD(node-1,kdof)=dval;
    CHKERRQ(1);
  }
  PetscPrintf(PETSCFEM_COMM_WORLD,"Total loads: %d\n",nloads);
#endif

  //fclose(fid);
  astr_destroy(linecopy);

  // This is never called!
  // It's just a trick to ensure that the
  // `parallel_read' routines are implemented in the binary
  static int dum=0;
  if (dum) {
#define DOINST(type)				\
    dvector<type> v_##type;			\
    dvector_read_parallel(NULL,v_##type)
    DOINST(int);
    DOINST(double);
    DOINST(char);
    DOINST(float);
  }
  return 0;
}

