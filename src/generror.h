// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: generror.h,v 1.5 2003/12/06 16:00:25 mstorti Exp $
#ifndef PETSCFEM_GENERROR_H
#define PETSCFEM_GENERROR_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class GenericError : public string { 
public:
  GenericError() : string("") { }
  GenericError(string s) : string(s) { }
  // GenericError(char *s) : string(s) { }
  GenericError(char *s,va_list l);
  GenericError(char *s,...);
};

#define PETSCFEM_ASSERT_GE(cond,templ,...)			\
if (!(cond)) { throw GenericError(templ "\n---------------\n"	\
	      "PETSC-FEM error at file \"%s\", line %d\n",	\
	       __VA_ARGS__,__FILE__,__LINE__); }

#define PETSCFEM_ASSERT_GE0(cond,templ) 				\
     if (!(cond)) { throw GenericError(templ "\n---------------\n"	\
     "PETSC-FEM error at file \"%s\", line %d\n",			\
     __FILE__,__LINE__); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void 
petscfem_check_par_err(int ierro, GenericError &e);

extern GenericError PETSCFEM_GENERIC_ERROR;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define CHECK_PAR_ERR_GE			\
catch(GenericError e) { ierro = 1; PETSCFEM_GENERIC_ERROR=e; }	\
petscfem_check_par_err(ierro,PETSCFEM_GENERIC_ERROR);

#endif
