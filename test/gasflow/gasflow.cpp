//__INSERT_LICENSE__
//$Id: gasflow.cpp,v 1.1 2003/01/25 21:25:49 mstorti Exp $
#define _GNU_SOURCE

extern int MY_RANK,SIZE;
#include <src/ampli.h>
#include <cstdio>
#include <sys/stat.h>
#include <cassert>
#include <cstdlib>
#include <vector>

#include <src/vecmacros.h>
#include <src/fstack.h>
#include <src/texthash.h>
#include <src/fem.h>
#include <src/texthf.h>

#include <src/hook.h>
#include <src/dlhook.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class gasflow_hook {
public:
  void init(Mesh &mesh,Dofmap &dofmap,
	    TextHashTableFilter *options,const char *name) {
    printf("gasflow_hook: init\n"); }
  void time_step_pre(double time,int step) {
    printf("gasflow_hook: time_step_pre\n"); }
  void time_step_post(double time,int step,
		      const vector<double> &gather_values) {
    printf("gasflow_hook: time_step_post\n"); }
  void close() {
    printf("gasflow_hook: close\n"); }
};

DL_GENERIC_HOOK(gasflow_hook);
