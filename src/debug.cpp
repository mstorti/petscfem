//__INSERT_LICENSE__
//$Id: debug.cpp,v 1.4 2001/11/22 19:50:11 mstorti Exp $
 
#include <cstdio>
#include <time.h>

#include <petsc.h>
#include <src/debug.h>

Debug debug;

int Debug::active(const char *s=NULL) const {
  string ss = string(s!=NULL ? s : "");
  map<string,int>::const_iterator k;
  k = active_flags.find(ss);
  if (k == active_flags.end() || ! k->second) {
    return 0;
  } else {
    return 1;
  }
}

void Debug::activate(const char *s=NULL) {
  string ss = string(s!=NULL ? s : "");
  active_flags[ss] = 1;
}

void Debug::deactivate(const char *s=NULL) {
  string ss = string(s!=NULL ? s : "");
  active_flags[ss] = 0;
}

void Debug::trace(const char *s=NULL) {
  time_t tt;
#define MXTM 100
  char t[MXTM];
  MPI_Comm_rank(comm,&myrank);
  if ((active() || active("print")) && myrank==0) {
    tt = time(NULL);
    // t = asctime(localtime(&tt));
    strftime(t,MXTM,"%H:%M:%S",localtime(&tt));
    printf("-- %s -- [%s %10.3f]",s,t,chrono.elapsed());
  }
  if (!active()) {
    printf("\n");
    return;
  }
  int ierr;
  char ans;
  ierr = MPI_Barrier(comm);
  assert(ierr==0);
  if (myrank==0) {
    printf("Continue? (n/RET=y) > ");
    fflush(stdout);
    scanf("%c",&ans);
  }
  ierr = MPI_Bcast (&ans, 1, MPI_CHAR, 0,comm);
  assert(ierr==0); 
  if (ans=='n') {
    PetscFinalize();
    exit(0);
  } else if (ans=='d') {
    deactivate();
  } 
  ierr = MPI_Barrier(comm);
  assert(ierr==0);  
}

Debug::Debug(int active_=0,MPI_Comm comm_=MPI_COMM_WORLD) : 
  comm(comm_) {
  if (active_) activate();
  chrono.start();
}
