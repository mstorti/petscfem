// -*- mode: c++ -*-
//__INSERT_LICENSE__
//$Id: genload.h,v 1.3 2001/04/14 13:21:03 mstorti Exp $

#ifndef GENLOAD_H
#define GENLOAD_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class genload : public Elemset {
public:
  int assemble(double *retval,NodeData *nodedata,double *locst,
	       double *locst2,Dofmap *dofmap,int ijob,
	       char *jobinfo,int myrank,int el_start,int el_last,int iter_mode);
};

#endif
