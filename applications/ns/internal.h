/*__INSERT_LICENSE__*/
//$Id: internal.h,v 1.5 2001/05/30 03:58:44 mstorti Exp $
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
 
