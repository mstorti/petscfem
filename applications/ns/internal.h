//__INSERT_LICENSE__
//$Id: internal.h,v 1.4 2001/05/02 00:08:59 mstorti Exp $
#ifndef INTERNAL_H
#define INTERNAL_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class internal : public Elemset {
public:
  int assemble(double *retval,Nodedata *nodedata,double *locst,
	       Dofmap *dofmap,int ijob,int myrank,int el_start,
	       int el_last,int iter_mode);
};

#endif

//Local Variables: *
//mode: c++ *
//End: *
 
