// -*- mode: c++ -*-
/*__INSERT_LICENSE__*/
//$Id: fem.h,v 1.39.10.1 2007/02/19 20:23:56 mstorti Exp $

#ifndef FEM_H
#define FEM_H

using namespace std;

#include <newmatio.h>
#include <petscksp.h>
#include <stdlib.h>

#if !defined(PETSC_VERSION_)
#define PETSC_VERSION_(MAJOR,MINOR,SUBMINOR) \
  ((PETSC_VERSION_MAJOR == (MAJOR)) &&       \
   (PETSC_VERSION_MINOR == (MINOR)) &&       \
   (PETSC_VERSION_SUBMINOR == (SUBMINOR)) && \
   (PETSC_VERSION_RELEASE  == 1))
#endif

extern MPI_Comm PETSCFEM_COMM_WORLD;
extern int MY_RANK,SIZE;
int PetscFemInitialize(int *argc,char ***args,const char file[],const char help[]);


// Libretto has some incompatibilities with recent versions
// of compilers, so I have a fixed version of the header in ./src
#include "libretto.h"
#include <libretto/darray.h>

#undef HAVE_MEMMOVE // para que no chille
#include <src/gpdata.h>
#include <src/texthash.h>
#include <src/vecmacros.h>
#include <src/dofmap.h>
#include <src/buffpack.h>
#include <src/generror.h>

/**@name fem.h */
//@{
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** If user enters an elemset type #alias# in the data
    file, associate with real elemset #elemset_type#. 
    @author M. Storti
    @param elem_set_type_alias elemset type alias
    @param elem_set_type real elemset type
*/ 
#define SET_ELEMSET_TYPE_ALIAS(alias,elemset_type)	\
      if ( ! strcmp(type,#alias)) {			\
  	elemset = (Elemset *)new elemset_type;		\
      } else 

/** Use this macro to include elemset types in the function
    bless\_elemset. 
    NOTE: This is somewhat incorrect. NewElemset in the future will
    will replace Elemset. An in that case we will no need the explicit
    caset #(Elemset *)#. This is due to the fact that
    #Elemset# has private members for NewElemset and derived
    classes, so that we can not make polymorphism between these
    classes. 
    @author M. Storti
    @param elemset_type new elemset type to be included
*/ 
#define SET_ELEMSET_TYPE(elemset_type)			\
        SET_ELEMSET_TYPE_ALIAS(elemset_type,elemset_type)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This is prepended in some of the GETOPT macros to the name of the
    variable. 
    Remember to first "undef" the macro
*/ 
#define GETOPTDEF_HOOK(name) name

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Gets a value of type int or double from any hash table. 
    Example: TGETOPTDEF(thash,int,n,10) 
    @author M. Storti
    @param thash the TextHashTable from where to get the value
    @param type may be `int' or `double'
    @param name name of the variable
    @param default default value. 
*/ 
#define TGETOPTDEF(thash,type,name,default)			\
        type name=default;					\
        ierr = get_##type(thash,#name,&name,1);			\
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Gets a value of type int or double from any hash table. 
    No default value assumed.
    Example: GETOPTDEF(thash,int,n,none) 
    @author M. Storti
    @param thash the TextHashTable from where to get the value
    @param type may be `int' or `double'
    @param name name of the variable
    @param default (none) for use with `odoc.pl' default value. 
*/ 
#define TGETOPTNDEF(thash,type,name,default)			\
        type name;						\
        ierr = get_##type(thash,#name,&name);			\
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Gets a value of type string. 
    Example: {\tt TGETOPTDEF\_S(thash,string,preco\_type,"Jacobi")}
    @author M. Storti
    @param thash the TextHashTable from where to get the value
    @param type should be string
    @param name name of the variable
    @param default default value (not in parentheses)
*/ 
#define TGETOPTDEF_S(thash,type,name,default)			\
        type name=type(#default);				\
        ierr = get_##type(thash,#name,name,1);			\
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Gets a value of type string (doesn't declare the variable). 
    Example: {\tt TGETOPTDEF\_S\_ND(thash,string,preco\_type,"Jacobi")}
    @author M. Storti
    @param thash the TextHashTable from where to get the value
    @param type should be string
    @param name name of the variable (should be declared previously)
    @param default default value (not in parentheses)
*/ 
#define TGETOPTDEF_S_ND(thash,type,name,default)		\
        name=type(#default);					\
        ierr = get_##type(thash,#name,name,1);			\
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Gets a value of type string from the thash of the elemset
    Example: {\tt NTGETOPTDEF\_S(string,preco\_type,"Jacobi")}
    @author M. Storti
    @param type should be string
    @param name name of the variable
    @param default default value (not in parentheses)
*/ 
#define NGETOPTDEF_S(type,name,default) \
        type name=type(#default); \
        ierr = get_##type(#name,name,1); \
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Gets a value of type int or double from the global hash table. 
    Example: GETOPTDEF(int,n,10) 
    @author M. Storti
    @param type may be `int' or `double'
    @param name name of the variable
    @param default default value. 
*/ 
#define GETOPTDEF(type,name,default)				\
        type name=default;					\
        ierr = get_##type(mesh->global_options,#name,&name,1);	\
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Gets a value of type int or double from the global options hash table. 
    Example: GGETOPTDEF(int,n,10) 
    @author M. Storti
    @param type may be `int' or `double'
    @param name name of the variable
    @param default default value. 
*/ 
#define GGETOPTDEF(type,name,default)				\
             TGETOPTDEF(GLOBAL_OPTIONS,type,name,default)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Gets a value of type int or double from the general
    element hash table. 
    Example: SGETOPTDEF(int,n,10)
    @author M. Storti
    @param type may be `int' or `double'
    @param name name of the variable
    @param default default value. 
*/ 
#define SGETOPTDEF(type,name,default)				\
        type name=default;					\
        ierr = get_##type(thash,#name,&name,1);			\
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Gets a value of type int or double from the elemset
    hash table. Doesn't define the variable.
    Example: SGETOPTDEF(int,n,10)
    @author M. Storti
    @param type may be `int' or `double'
    @param name name of the variable
    @param default default value. 
*/ 
#define SGETOPTDEF_ND(type,name,default)				\
        (GETOPTDEF_HOOK(name)) = default;				\
        ierr = get_##type(thash,#name,&(GETOPTDEF_HOOK(name)),1);	\
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Gets a value of type int or double from an elemset
    hash table. 
    Example: EGETOPTDEF(elemset,int,n,10)
    @author M. Storti
    @param type may be `int', `double' or 'string'
    @param name name of the variable
    @param default default value. 
*/ 
#define EGETOPTDEF(elemset,type,name,default)			\
        type name;						\
        name = default;						\
        ierr = elemset->get_##type(#name,name,1);		\
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Gets a value of type int or double from an elemset
    hash table. Doesn't define the variable.
    Example: EGETOPTDEF_ND(elemset,int,n,10)
    @author M. Storti
    @param type may be `int', `double' or 'string'
    @param name name of the variable
    @param default default value. 
*/ 
#define EGETOPTDEF_ND(elemset,type,name,default)		\
        name = default;						\
        ierr = elemset->get_##type(#name,name,1);		\
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Gets a value of type int or double from any 
    hash table. Doesn't define the variable.
    Example: GGETOPTDEF(thash,int,n,10)
    @author M. Storti
    @param thash the TextHashTable from where to get the value
    @param type may be `int' or `double'
    @param name name of the variable
    @param default default value. 
*/ 
#define TGETOPTDEF_ND(thash,type,name,default) \
        name = default; \
        ierr = get_##type(thash,#name,&name,1); \
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Gets a value of type int or double from the general
    element hash table. OOP version. 
    Example: NSGETOPTDEF(int,n,10)
    @author M. Storti
    @param type may be `int' or `double'
    @param name name of the variable
    @param default default value. 
*/ 
#define NSGETOPTDEF(type,name,default)				\
        type name=default;					\
        ierr = get_##type(#name,name,1);			\
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Gets a value of type int or double from the general
    element hash table. OOP version. Doesn't define the variable.
    Example: NSGETOPTDEF_ND(int,n,10)
    @author M. Storti
    @param type may be `int' or `double'
    @param name name of the variable
    @param default default value. 
*/ 
#define NSGETOPTDEF_ND(type,name,default)			\
        GETOPTDEF_HOOK(name)=default;				\
        ierr = get_##type(#name,GETOPTDEF_HOOK(name),1);	\
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Sets an error condition and back-traces (routines).
    Based on PETSc `PFEMERRQ' macro. Use this macro in internal
    routines. Use PFEMERRA in the main.  
    @author M. Storti
    @param s string of error message
*/ 
#define PFEMERRQ(s) {PetscPrintf(PETSCFEM_COMM_WORLD,s); CHKERRQ(1);}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Sets an error condition and back-traces (main).
    Based on PETSc `PFEMERRA' macro. Use this macro in the main.
    Use PFEMERRQ in internal routines. 
    @author M. Storti
    @param s string of error message
*/ 
#define PFEMERRA(s) {PetscPrintf(PETSCFEM_COMM_WORLD,s); CHKERRA(1);}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Sets an error condition depending on error code and back-traces
    (routines).  Like PFEMERRQ, but checks integer `ierr' and sets
    error if ierr>0.  Based on PETSc `PFEMERRQ' macro. Use this macro
    in internal routines. Use PFEMERRA in the main.

    @author M. Storti
    @param ierr error code (input)
    @param s string of error message */
#define PFEMERRCQ(ierr,s) if (ierr)				\
         {PetscPrintf(PETSCFEM_COMM_WORLD,s); CHKERRQ(ierr);}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Sets an error condition depending on error code and back-traces
    (routines).  Like PFEMERRQ, but checks integer `ierr' and sets
    error if ierr>0.  Based on PETSc `PFEMERRQ' macro. Use this macro
    in the main. Use PFEMERRA in internal routines.

    @author M. Storti
    @param ierr error code (input)
    @param s string of error message */
#define PFEMERRCA(ierr,s) if (ierr) {PetscPrintf(PETSCFEM_COMM_WORLD,s); assert(0);}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Issues an error (with #PetscPrintf(...))# and
    calls #PetscFinalize()#. Argument are as for 'printf()'
    usage: PETSCFEM_ERROR(template,args).
*/ 
#define USE_VARARG_MACROS

#ifdef USE_VARARG_MACROS
#define PETSCFEM_ERROR0(templ)                                  \
  petscfem_error(templ "\n---------------\n"                    \
                 "PETSC-FEM error at file \"%s\", line %d\n",	\
		 __FILE__,__LINE__)
#define PETSCFEM_ERROR(templ,...)				\
  petscfem_error(templ "\n---------------\n"			\
	      "PETSC-FEM error at file \"%s\", line %d\n",	\
		 __VA_ARGS__,__FILE__,__LINE__)

#define PETSCFEM_WARN0(templ)                                  \
  petscfem_warn(templ "\n---------------\n"                    \
                 "PETSC-FEM warning at file \"%s\", line %d\n",	\
		 __FILE__,__LINE__)
#define PETSCFEM_WARN(templ,...)				\
  petscfem_warn(templ "\n---------------\n"			\
	      "PETSC-FEM warning at file \"%s\", line %d\n",	\
		 __VA_ARGS__,__FILE__,__LINE__)
#else
#define PETSCFEM_ERROR petscfem_error
#define PETSCFEM_ERROR0 petscfem_error
#endif

void petscfem_printf(const char *templ,va_list list);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** If #bool_cond# evaluates to false issues 
    an error (with #PETSCFEM_ERROR#).
    usage: #PETSCFEM_ASSERT(bool_cond,printf_args)#;
*/ 
#ifdef USE_VARARG_MACROS
// Bool condition and message
#define PETSCFEM_ASSERT0(bool_cond,templ)       \
if (!(bool_cond)) {                             \
  PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,  \
                  "Assertion failed: \"%s\"\n", \
                   #bool_cond);                 \
  PetscSynchronizedFlush(PETSCFEM_COMM_WORLD);  \
  PETSCFEM_ERROR0(templ);}

// Only bool condition
#define PETSCFEM_ASSERTB(bool_cond) \
  PETSCFEM_ASSERT0(bool_cond,"")

// Bool condition and message with args 
#define PETSCFEM_ASSERT(bool_cond,templ,...)    \
if (!(bool_cond)) {                             \
  PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,  \
                  "Assertion failed: \"%s\"\n", \
                  #bool_cond);                  \
  PetscSynchronizedFlush(PETSCFEM_COMM_WORLD);  \
  PETSCFEM_ERROR(templ,__VA_ARGS__);}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// WARNING MACROS
// Only bool condition
#define PETSCFEM_WARNINGB(bool_cond)            \
if (!(bool_cond)) {                             \
  printf("Warning (failed assertion): ");         \
  PETSCFEM_WARN0(#bool_cond);}

// Bool condition and message
#define PETSCFEM_WARNING0(bool_cond,templ)      \
if (!(bool_cond)) {                             \
  PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,  \
    "Warning (failed assertion): \"%s\"\n",     \
    #bool_cond);                                \
  PetscSynchronizedFlush(PETSCFEM_COMM_WORLD);  \
  PETSCFEM_WARN0(templ);}

// Bool condition and message with args 
#define PETSCFEM_WARNING(bool_cond,templ,...)   \
if (!(bool_cond)) {                             \
  PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,  \
       "Warning (failed assertion): \"%s\"\n",  \
       #bool_cond);                             \
  PetscSynchronizedFlush(PETSCFEM_COMM_WORLD);  \
  PETSCFEM_WARN(templ,__VA_ARGS__);}
#else
#define PETSCFEM_ASSERT petscfem_assert
#define PETSCFEM_ASSERT0 petscfem_assert
#endif

void petscfem_warn(const char *templ,...);
void petscfem_error(const char *templ,...);
void petscfem_assert(int bool_cond,const char *templ,...);

#define PFEM_TRACE(s) PetscPrintf(PETSCFEM_COMM_WORLD,		\
     "<%s>. At file " __FILE__ ", line %d\n",s,__LINE__)

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/** Contains (constant) data relative to nodes (may be coordinates and
    other). 
    @author M. Storti
    @param nnod number of nodes
    @param nu number of real quantities per node
*/ 
class Nodedata {
public:
  ~Nodedata();
  Nodedata();
  double *nodedata;
  int nnod,ndim,nu;
  typedef pair<int,double*> Field;
  typedef map<string,Field> FieldMap;
  FieldMap fields;
  TextHashTable *options;
  bool get_field(const string&,int*,double**);
};

class Elemset;
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Contains a list of the elemsets, node information and a hash table
    of global options. 
    @author M. Storti
    @param elemsetlist list of elemsets
    @param nodedata constant data per node
    @param global\_options hash table with global options
*/ 
class Mesh {
public:
  /// The list of elemsets
  Darray *elemsetlist;
  /// The array of nodes
  Nodedata *nodedata;
  /// The option table
  TextHashTable *global_options;
  /// Finds an elemset give its name
  Elemset *find(const string &name);
  /// Default ctor
  Mesh();
  // dtor
  ~Mesh();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

#include "elemset.h"

// Other defs
#define FLEN 300  // length of filename
#define ELEM_CHUNK_SIZE 1000

// Values for ijob
/** @name ijob
    Determines what kind of action is performed on the data
*/ 
#define COMP_MAT 1  // Compute matrices
#define COMP_VEC 2  // Compute vector
#define COMP_FDJ 3  // Compute finite difference jacobian
#define COMP_FDJ_PROF 4  // Compute finite difference jacobian profile
#define COMP_MAT_PROF 5  // Compute finite difference jacobian profile

/// parametros numericos
#define EPSILON_FDJ 1e-6
//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/** @name Utility definitions */
//@{
/// Prints a variable along with his name (useful for debugging). 
#define SHV_OPTIONS
#define SHV(x) cout << #x ": " SHV_OPTIONS << x << endl

/// Prints a matrix variable along with his name (useful for debugging). 
#define SHM(x) cout << #x ": " << x << endl

/// Prints a trace
#define TRACE(s) PetscPrintf(PETSCFEM_COMM_WORLD,">>>TRACE <" \
                   #s ">  (%s:%d)\n",__FILE__,__LINE__)

/// Prints a variable with printf()
#define SHVS(name,f)  printf(#name ": %" #f "\n",name)
//@}
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
//int get_line(FILE *fid,char *line);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Converts a state vector (reduced form) to node/field form and
    prints it. 
    @author M. Storti
    @param filename (input) file where to write the vector. May contain
    relative directories. 
    @param x (input) PETSc MPI vector to be written 
    @param dofmap (input) corresponding dofmap 
    @param time_data (input, def=NULL) an external parameter in order to compute
    external boundary conditions, etc...
    @param append (input, def=0) appending mode (append if
    `append==0') 
*/ 
int print_vector(const char *filename,const Vec x,
		 const Dofmap *dofmap,const TimeData *time_data=NULL,
		 const int append=0);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Converts a state vector (reduced form) to node/field form. 
    @author M. Storti
    @param fields (output) double array where the node/fields form is put. 
    @param x (input) PETSc MPI vector to be written 
    @param dofmap (input) corresponding dofmap 
    @param time_data (input, def=NULL) an external parameter in order to compute
    external boundary conditions, etc...
    @return PETSc error code
*/ 
int state2fields(double *fields,const Vec x,const Dofmap *dofmap,
		 const TimeData *time_data=NULL);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Reads a vector from a file.
    This can be used for initialization for instance. 
    @author M. Storti
    @param filename file from where to read  the vector. May contain
    relative directories. 
    @param x PETSc MPI vector to be read
    @param x PETSc sequantial working vector
    @param dofmap corresponding dofmap 
*/ 
int read_vector(const char *filename,Vec x,Dofmap *dofmap,int
		myrank);

// obsolete (???)
// int read_vector2(char *filename,Vec x,Vec xseq,Dofmap *dofmap,int myrank);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Computes the sparse profile for defining matrices.  PETSc needs
    the number of non null entries in the diagonal and off diagonal
    blocks. Calling `assemble' with option `COMP\_MAT\_PROF' computes
    `da' and this function defines a matrix prototype. 

    @author M. Storti
    @param da (input) the dynamic array computed by assemble
    @param dofmap (input) the corresponding dofmap
    @param myrank (input) the index of the current processor
    @param A (output) PETSc matrix prototype
*/
int compute_prof(Darray *da,Dofmap *dofmap,int myrank, Mat *A,
		 int debug_compute_prof);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
/** Sets a Matrix to zero.
    This is obsolete. We should call MatZeroEntries() (PETSc)

    @author M. Storti
    @param A matrix to set to zero
    @param ass_flag a flag indicating wether the matrix was already
    built or not. 
*/ 
int zeroe_mat(Mat A,int & ass_flag);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** @memo Optionally reads a vector if `initial\_name <string>' is in the
    global\_options hash text table. If not, sets initial vector to 0. 
    @author M. Storti
    @param mesh (input) the mesh read. 
    @param x (input) the vector to be read. 
    @param dofmap (input) the dofmap of the mesh. 
    @param myrank (input) the number of this processor. 
*/ 
int opt_read_vector(Mesh *mesh,Vec x, Dofmap *dofmap,int myrank);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** @memo Prints copyright info. To be called at start time in the main
    program. 
    @author M. Storti
*/ 
void print_copyright(void);

#define GET_JOBINFO_FLAG(name) \
   int name  = !strcmp(jobinfo,#name)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Process this jobinfo for this elemset. 
    @author M. Storti
    @param name (input) the jobinfo to be processed
*/ 
#define DONT_SKIP_JOBINFO(name) \
   if (!strcmp(jobinfo,#name)) skip_elemset = 0

#define WAIT(s) wait_from_console(s "  --- at file: " __FILE__)

/** Returns the square of a double
    @param x (input) the base 
*/ 
inline double square(double x) { return x*x; }

/** Returns the cube of a double
    @param x (input) the base 
*/ 
inline double cube(double x) { return x*x*x; }

//#define SQ(n) ((n)*(n))
//#define CB(n) ((n)*(n)*(n))
#define SQ(n) square(n)
#define CB(n) cube(n)

// Each module has a different of this that prints the
// time at which it was compiled, compiler version, etc... 
// Normally contained in a file named `compadte.cpp'. 
void print_petscfem_link_date(void);

// Prints the day of the run. 
void print_date(void);

// Signals an error with computation of Jacobian
void detj_error(double &detJaco,int elem);

// CHKERRA is now obsolete in PETSc (see changes for PETSc 2.1)
#define CHKERRA CHKERRQ

#if PETSC_VERSION_(3,1,0) || PETSC_VERSION_(3,0,0)

typedef PetscTruth PetscBool;

#define PETSCVIEWERASCII PETSC_VIEWER_ASCII

#define ISCreateGeneral(a,b,c,d,e)\
        ISCreateGeneral(a,b,c,e)

inline PetscErrorCode PetscViewerDestroy_Compat(PetscViewer *a)
{PetscViewer b = *a; *a=0; return PetscViewerDestroy(b);}
#define PetscViewerDestroy PetscViewerDestroy_Compat
inline PetscErrorCode ISDestroy_Compat(IS *a)
{IS b = *a; *a=0; return ISDestroy(b);}
#define ISDestroy ISDestroy_Compat
inline PetscErrorCode VecScatterDestroy_Compat(VecScatter *a)
{VecScatter b = *a; *a=0; return VecScatterDestroy(b);}
#define VecScatterDestroy VecScatterDestroy_Compat
inline PetscErrorCode VecDestroy_Compat(Vec *a)
{Vec b = *a; *a=0; return VecDestroy(b);}
#define VecDestroy VecDestroy_Compat
inline PetscErrorCode MatDestroy_Compat(Mat *a)
{Mat b = *a; *a=0; return MatDestroy(b);}
#define MatDestroy MatDestroy_Compat
inline PetscErrorCode KSPDestroy_Compat(KSP *a)
{KSP b = *a; *a=0; return KSPDestroy(b);}
#define KSPDestroy KSPDestroy_Compat

#define KSPSetPCSide KSPSetPreconditionerSide

#endif


#if PETSC_VERSION_(2,3,2)
#define VecScatterBegin(ctx,x,y,im,sm) \
        VecScatterBegin((x),(y),(im),(sm),(ctx))
#define VecScatterEnd(ctx,x,y,im,sm) \
        VecScatterEnd((x),(y),(im),(sm),(ctx))
#endif

#if PETSC_VERSION_(2,3,2)
#define KSPMonitorSet KSPSetMonitor
#define KSPMonitorCancel KSPClearMonitor
#define KSPMonitorDefault  KSPDefaultMonitor
#define KSPMonitorTrueResidualNorm KSPTrueMonitor
#define KSPMonitorSolution KSPVecViewMonitor
#define KSPMonitorLG KSPLGMonitor
#endif

#if PETSC_VERSION_(2,3,3) || PETSC_VERSION_(2,3,2)
#define MatSetOption(mat,opt,flg) MatSetOption( (mat), (assert(flg),(opt)) )
#endif

#define PFUNUSED __attribute__ ((unused)) 

#endif
//@}
