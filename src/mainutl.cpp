//__INSERT_LICENSE__
//$Id: mainutl.cpp,v 1.22 2005/02/20 16:55:37 mstorti Exp $
 
#include "fem.h"
#include "utils.h"
#include "util2.h"
#include "readmesh.h"
#include "idmap.h"
#include "elemset.h"
#include "generror.h"
#include <src/dvector.h>

extern int MY_RANK,SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "readval(int&,char *,double&)"
int readval(int &rflag,char *line,double &val) {
  static char *bsp=" \t";
  char *token;
  token = strtok((rflag==0 ? rflag=1,line : NULL),bsp);
  if (token==NULL) return 1;
  int nread = sscanf(token,"%lf",&val);
  return nread!=1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "readval(int&,char *,int&)"
int readval(int &rflag,char *line,int &val) {
  static char *bsp=" \t";
  char *token;
  token = strtok((rflag==0 ? rflag=1,line : NULL),bsp);
  int nread = sscanf(token,"%d",&val);
  return nread!=1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "read_hash_table"
#if 1
int read_hash_table(FileStack *& fstack,TextHashTable *& thash) {
  thash = new TextHashTable;
  thash->read(fstack);
  return 0;
}
#else
int read_hash_table(FileStack *& fstack,TextHashTable *& thash) {

  char *line;
  char *key, *val,*bsp=" \t";
  thash = new TextHashTable;
  int he=0;
  while (1) {
    fstack->get_line(line);
    if (strstr("__END_HASH__",line)) break;
    key = strtok(line,bsp);
    val = strtok(NULL,"\n");
    if (!strcmp(key,"_table_include")) {
      thash->include_table(val);
    } else {
      he++;
      if (val==NULL) val="";
      thash->add_entry(key,val);
    }
  }
  g_hash_table_freeze(thash->hash);
  return 0;
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "print_vector_rota" 
void print_vector_rota(const char *filenamepat,const Vec x,const
		       Dofmap *dofmap,const TimeData *time_data,
		       const int j,const int nsave,const int nrec,
		       const int nfile) {
  
  if (nsave==0) return;
  div_t ressave,resrec;
  ressave = div(j,nsave);
  if (ressave.rem !=0) return;

  char buf[300];
  resrec = div(ressave.quot,nrec);
  int irec=resrec.rem;
  int ifile = resrec.quot % nfile;

  sprintf(buf,filenamepat,ifile);
  PetscPrintf(PETSC_COMM_WORLD,
	      "print_vector_rota: step = %d, saving on rec %d, file %s\n",j,irec,buf);
  print_vector(buf,x,dofmap,time_data,irec!=0);

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "print_vector" 
int print_vector(const char *filename,const Vec x,const Dofmap *dofmap,
		 const TimeData *time_data,const int append) {

  double *vseq_vals,*sstate;
  Vec vseq;
  
  int myrank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);

  // fixme:= Now we can make this without a scatter. We can use
  // the version of get_nodal_value() with ghost_values. 
  int neql = (myrank==0 ? dofmap->neq : 0);
  int ierr = VecCreateSeq(PETSC_COMM_SELF,neql,&vseq);  CHKERRQ(ierr);
  ierr = VecScatterBegin(x,vseq,INSERT_VALUES,
			 SCATTER_FORWARD,*dofmap->scatter_print); CHKERRA(ierr); 
  ierr = VecScatterEnd(x,vseq,INSERT_VALUES,
		       SCATTER_FORWARD,*dofmap->scatter_print); CHKERRA(ierr); 
  ierr = VecGetArray(vseq,&vseq_vals); CHKERRQ(ierr);
 
  if (myrank==0) try {
    printf("Writing vector to file \"%s\"\n",filename);
    FILE *output;
    output = fopen(filename,(append == 0 ? "w" : "a" ) );
    if (output==NULL) throw GenericError("Couldn't open output file");

    int ndof=dofmap->ndof;
    double dval;
    for (int k=1; k<=dofmap->nnod; k++) {
      for (int kldof=1; kldof<=ndof; kldof++) {
	dofmap->get_nodal_value(k,kldof,vseq_vals,time_data,dval);
	fprintf(output,"%12.10e  ",dval);
      }
      fprintf(output,"\n");
    }
    fclose(output);
  } catch(GenericError e) { ierr = 1; }
  MPI_Bcast(&ierr,1,MPI_INT,0,PETSC_COMM_WORLD);
  PETSCFEM_ASSERT0(!ierr,"Couldn't open output file");
  ierr = VecRestoreArray(vseq,&vseq_vals); CHKERRQ(ierr); 
  ierr = VecDestroy(vseq);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "state2fields" 
int state2fields(double *fields,const Vec x,const Dofmap *dofmap,
		 const TimeData *time_data) {

  double *vseq_vals,*sstate;
  Vec vseq;
  
  // fixme:= Now we can make this without a scatter. We can use
  // the version of get_nodal_value() with ghost_values. 
  int neql = (!MY_RANK ? dofmap->neq : 0);
  int ierr = VecCreateSeq(PETSC_COMM_SELF,neql,&vseq);  CHKERRQ(ierr);
  ierr = VecScatterBegin(x,vseq,INSERT_VALUES,
			 SCATTER_FORWARD,*dofmap->scatter_print); CHKERRA(ierr); 
  ierr = VecScatterEnd(x,vseq,INSERT_VALUES,
		       SCATTER_FORWARD,*dofmap->scatter_print); CHKERRA(ierr); 
  ierr = VecGetArray(vseq,&vseq_vals); CHKERRQ(ierr);
 
  if (!MY_RANK) {
    int ndof=dofmap->ndof;
    double dval;
    for (int k=1; k<=dofmap->nnod; k++) {
      for (int kldof=1; kldof<=ndof; kldof++) {
	dofmap->get_nodal_value(k,kldof,vseq_vals,time_data,dval);
	fields[ndof*(k-1)+kldof-1] = dval;
      }
    }
  }
  ierr = VecRestoreArray(vseq,&vseq_vals); CHKERRQ(ierr); 
  ierr = VecDestroy(vseq);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "print_some(const char *,const State &,Dofmap *,set<int> &)" 
int print_some(const char *filename,const State &s,Dofmap *dofmap,
	       set<int> & node_list) {
  print_some(filename,s.v(),dofmap,node_list,&s.t());
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "print_some(const char *,const Vec ,Dofmap *," \
"set<int> ,const TimeData *=NULL)" 
int print_some(const char *filename,const Vec x,Dofmap *dofmap,
	       set<int> node_list,const TimeData *time_data) {

  Vec vseq;
  double *sol;
  int ierr;

  int myrank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);

  // fixme:= Now we can make this without a scatter. We can use
  // the version of get_nodal_value() with ghost_values. 
  int neql = (myrank==0 ? dofmap->neq : 0);
  ierr = VecCreateSeq(PETSC_COMM_SELF,neql,&vseq);  CHKERRQ(ierr);

  ierr = VecScatterBegin(x,vseq,INSERT_VALUES,
			 SCATTER_FORWARD,
			 *(dofmap->scatter_print)); CHKERRA(ierr); 
  ierr = VecScatterEnd(x,vseq,INSERT_VALUES,
		       SCATTER_FORWARD,
		       *(dofmap->scatter_print)); CHKERRA(ierr); 

  ierr = VecGetArray(vseq,&sol); CHKERRA(ierr); 
  if (myrank==0) {
    printf("Printing some to file \"%s\"\n",filename);
    FILE *output;
    output = fopen(filename,"a");
    if (output==NULL) {
      printf("Couldn't open output file\n");
      // fixme:= esto esta mal. Todos los procesadores
      // tienen que llamar a PetscFinalize()
      PetscFinalize();
      exit(1);
    }

    int ndof=dofmap->ndof;
    double dval;
    for (int k=1; k<=dofmap->nnod; k++) {
      if (node_list.find(k)!=node_list.end()) {
	fprintf(output,"%d   ",k);
	for (int kldof=1; kldof<=ndof; kldof++) {
	  dofmap->get_nodal_value(k,kldof,sol,time_data,dval);
	  fprintf(output,"%12.10e  ",dval);
	}
	fprintf(output,"\n");
      }
    }
    fclose(output);
  }
  ierr = VecRestoreArray(vseq,&sol); CHKERRA(ierr); 
  ierr = VecDestroy(vseq);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "read_vector" 
int read_vector(const char *filename,Vec x,Dofmap *dofmap,int myrank) {

  int ndof = dofmap->ndof;
  int ierr,code,warn_flag=0,ierro=0;

  PetscPrintf(PETSC_COMM_WORLD,"Reading vector from file \"%s\"\n",filename);
  dvector<double> xdof(dofmap->neqtot);
  xdof.resize(dofmap->neqtot);
  if (myrank==0) {
    dvector<double> xext(dofmap->nnod*ndof);
    xext.a_resize(2,dofmap->nnod,ndof);
    FILE *fid;
    fid = fopen(filename,"r");
    if (fid==NULL) {
      printf("read_vector: couldn't open file %s\n",
	     filename);
      ierro=1;
      throw ("Cant open file");
    }
    double dval;
    for (int k=1; k<=dofmap->nnod; k++) {
      // if (k % 1000 == 0) printf("%d node values read...\n",k);
      for (int kldof=1; kldof<=ndof; kldof++) {
	code = fscanf(fid,"%lf",&xext.e(k-1,kldof-1));
	if (code==EOF) {
	  warn_flag=1;
	  xext.e(k-1,kldof-1) = 0.;
	}
      }
    }
    fclose(fid);
    // printf("%d node values read...\n",dofmap->nnod);
    if (warn_flag) 
      PetscPrintf(PETSC_COMM_WORLD,
		  "PETScFEM warning: not enough values"
		  " read while reading vector. Filling with 0's.\n");
    dofmap->solve(xdof.buff(),xext.buff());
    xext.clear();
  } 
  CHECK_PAR_ERR(ierro,"Error reading vector from file.");

  ierr = MPI_Bcast (xdof.buff(),dofmap->neqtot,
		    MPI_DOUBLE,0,PETSC_COMM_WORLD);

  for (int k=0; k<dofmap->neq; k++) {
    if (dofmap->dof1 <= k+1 <= dofmap->dof2) {
      VecSetValue(x,k,xdof.e(k),INSERT_VALUES);
    }
  }
  xdof.clear();
  // PetscPrintf(PETSC_COMM_WORLD,"Values set.\n",dofmap->nnod);
  ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Done.\n",filename);
  return ierr;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int print_some_file_init(TextHashTable *thash," \
"const char *,const char *,set<int> &)"
int print_some_file_init(TextHashTable *thash,
			 const char *print_some_file,
			 const char *save_file_some,set<int> &node_list,
			 int save_file_some_append) {
  if (MY_RANK==0 && strlen(print_some_file)>0) {
    int nodo;
    PetscPrintf(PETSC_COMM_WORLD,"Reading print_some_file...\n");
    FILE *fid=fopen(print_some_file,"r");
    if (fid==NULL) {
      PetscPrintf(PETSC_COMM_WORLD,"Couldn't open `print_some_file': \"%s\"\n",
		  print_some_file);
      PetscFinalize();
      exit(0);
    }
    while (1) {
      int nread = fscanf(fid,"%d",&nodo);
      if (nread==EOF) break;
      node_list.insert(nodo);
    }
    fclose(fid);
    PetscPrintf(PETSC_COMM_WORLD,"... Done.\n");
    if (!save_file_some_append) {
      // Rewind file, discard old content
      FILE *fid = fopen(save_file_some,"w");
      fclose(fid);
    }
  }
  return 0;
}

