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
 
