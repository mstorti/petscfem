// -*- mode: c++ -*-
/*__INSERT_LICENSE__*/
//$Id: genload.h,v 1.5 2001/05/30 03:58:41 mstorti Exp $

#ifndef GENLOAD_H
#define GENLOAD_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class genload : public Elemset {
public:
  int assemble(double *retval,Nodedata *nodedata,double *locst,
	       double *locst2,Dofmap *dofmap,int ijob,
	       char *jobinfo,int myrank,int el_start,int el_last,int iter_mode);
};

#endif
