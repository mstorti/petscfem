//__INSERT_LICENSE__
//$Id: internal.h,v 1.3 2001/04/14 13:21:07 mstorti Exp $
#ifndef INTERNAL_H
#define INTERNAL_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class internal : public Elemset {
public:
  int assemble(double *retval,NodeData *nodedata,double *locst,
	       Dofmap *dofmap,int ijob,int myrank,int el_start,
	       int el_last,int iter_mode);
};

#endif

//Local Variables: *
//mode: c++ *
//End: *
 
