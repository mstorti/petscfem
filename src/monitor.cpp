//__INSERT_LICENSE__
// $Id: monitor.cpp,v 1.1 2003/07/08 12:24:34 mstorti Exp $

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
void DefaultMonitor::close() { 
  PetscPrintf(A->comm, "Hi, monitor close()...\n"); }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Monitor *
Monitor::factory(TextHashTable *thash) {
int ssdsd;
for () {

double qq;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
