//__INSERT_LICENSE__
//$Id: debug.cpp,v 1.1 2001/11/21 17:53:13 mstorti Exp $
 
#include <src/debug.h>

int Debug::active(const char *s=NULL) const {
  string ss = string(s!=NULL ? s : "");
  map<string,int>::iterator k;
  k = active_flags[ss];
  if (k != active_flags.end() || ! k->second) {
    return 0;
  } else {
    return 1;
  }
}

void activate(const char *=NULL) {
  string ss = string(s!=NULL ? s : "");
  active_flags[ss] = 1;
}

void activate(const char *=NULL) {
  string ss = string(s!=NULL ? s : "");
  active_flags[ss] = 0;
}

void trace(const char *s=NULL) {
  if (active() || active("print")) printf("%s --- ",s);
  if (!active()) return;
  MPI_Comm_rank(comm,&myrank);
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

Debug(int active_=0,MPI_Comm comm_=MPI_COMM_WORLD) : 
  comm(comm_) {
  if (active_) activate();
}


