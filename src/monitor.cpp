//__INSERT_LICENSE__
// $Id: monitor.cpp,v 1.4 2003/07/08 22:35:52 mstorti Exp $

#include <cstdio>
#include <petsc.h>
#include <src/monitor.h>
#include <src/monitor2.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DefaultMonitor::init(MPI_Comm comm_a,TextHashTable *options_a) {
  comm = comm_a;
  options = options_a;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DefaultMonitor::start() { 
  TGETOPTDEF(options,int,print_internal_loop_conv,1);
  PetscPrintf(comm, "Hi, monitor start()...\n"); 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DefaultMonitor::step(int n,double rnorm) { 
  if (print_internal_loop_conv) {
    if (n==0) PetscPrintf(comm,
			  " Begin internal iterations "
			  "--------------------------------------\n");
    PetscPrintf(comm,"iteration %d KSP "
		"Residual_norm = %14.12e \n",n,rnorm);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DefaultMonitor::stop() { 
  PetscPrintf(comm,"Hi, monitor close()...\n");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Monitor *Monitor::factory(TextHashTable *thash) {
  TGETOPTDEF_S(thash,string,monitor_type,default);
  Monitor * monitor=NULL;
//    if (monitor_type=="default") {
//      monitor = new DefaultMonitor;
  //    } else
  monitor->init(thash);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
