//__INSERT_LICENSE__
// $Id: monitor.cpp,v 1.2 2003/07/08 12:34:31 mstorti Exp $

#include <cstdio>
#include <petsc.h>
#include <src/monitor.h>
#include <src/monitor2.h>

void DefaultMonitor::start() { 
  PetscPrintf(A->comm, "Hi, monitor start()...\n"); 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DefaultMonitor::step(int n,double rnorm) { 
  if (A->print_internal_loop_conv) {
    if (n==0) PetscPrintf(A->comm,
			  " Begin internal iterations "
			  "--------------------------------------\n");
    PetscPrintf(A->comm,
		"iteration %d KSP Residual_norm = %14.12e \n",n,rnorm);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DefaultMonitor::stop() { 
  PetscPrintf(A->comm, "Hi, monitor close()...\n"); }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Monitor *Monitor::factory(TextHashTable *thash) {
TGETOPTDEF_S(thash,string,monitor_type,default);
Monitor * monitor=NULL;
if (monitor_type=="default") {
monitor = new DefaultMonitor;
} else {
monitor->init(thash);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
