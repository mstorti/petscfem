//__INSERT_LICENSE__
//$Id: fracstep.h,v 1.3 2001/04/14 13:21:07 mstorti Exp $
#ifndef FRACSTEP_H
#define FRACSTEP_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class fracstep : public Elemset {
public:
  int assemble(double *retval,NodeData *nodedata,double *locst,
	       double *locst2,Dofmap *dofmap,int ijob,
	       char *jobinfo,int myrank,int el_start,int el_last,int iter_mode);
};

#endif

//Local Variables: *
//mode: c++ *
//End: *
 
