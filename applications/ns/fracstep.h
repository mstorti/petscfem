//__INSERT_LICENSE__
//$Id: fracstep.h,v 1.2 2001/04/01 01:34:59 mstorti Exp $
#ifndef FRACSTEP_H
#define FRACSTEP_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class fracstep : public Elemset {
public:
  int assemble(double *retval,Nodedata *nodedata,double *locst,
	       double *locst2,Dofmap *dofmap,int ijob,
	       char *jobinfo,int myrank,int el_start,int el_last,int iter_mode);
};

#endif

//Local Variables: *
//mode: c++ *
//End: *
 
