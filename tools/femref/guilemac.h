// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: guilemac.h,v 1.1 2005/01/17 23:50:54 mstorti Exp $
#ifndef PETSCFEM_GUILEMAC_H
#define PETSCFEM_GUILEMAC_H

typedef SCM(*scm_fun)();

#define S(name) s_##name

#define MY_SCM_GET_ARG(name,tag,ctype,pos)	\
  SCM_ASSERT (SCM_SMOB_PREDICATE(tag,s_##name),	\
              s_##name, pos, FUNC_NAME);	\
  ctype name =					\
          (ctype)SCM_SMOB_DATA (s_##name)

#define DVDBLARG(name,pos)					\
  SCM_ASSERT(SCM_SMOB_PREDICATE(dvdbl_tag,s_##name),		\
              s_##name, pos, FUNC_NAME);			\
  dvector<double> *name =					\
          (dvector<double> *)SCM_SMOB_DATA (s_##name)

#define DVINTARG(name,pos)				\
   SCM_ASSERT (SCM_SMOB_PREDICATE(dvint_tag,s_##name),	\
              s_##name, pos, FUNC_NAME);		\
  dvector<int> *name =					\
          (dvector<int> *)SCM_SMOB_DATA (s_##name); 

#define MY_SCM_GET_INT_DEF(name,def,pos)			\
  int name;							\
  if (s_##name == SCM_UNDEFINED) name = def;			\
  else {							\
    SCM_ASSERT(SCM_INUMP(s_##name),s_##name,pos, __FUN__);	\
    name = SCM_INUM(s_##name);					\
  }

#define MY_SCM_GET_INT(name,pos)		\
  int name;					\
  SCM_ASSERT(SCM_INUMP(s_##name),		\
	     s_##name,pos, __FUN__);		\
  name = SCM_INUM(s_##name)

#define MY_SCM_GET_DBL(name,pos,msg)		\
  double name;					\
  SCM_ASSERT(scm_real_p(s_##name),		\
	     s_##name,pos, __FUN__);		\
  name = scm_num2dbl(s_##name,msg)

#define MY_SCM_GET_BOOL(name,pos)		\
  int name;					\
  SCM_ASSERT(scm_boolean_p(s_##name),		\
	     s_##name,pos, __FUN__);		\
  name = SCM_NFALSEP(s_##name)

#endif
