//__INSERT_LICENSE__
//$Id: debug.cpp,v 1.5.4.1 2001/12/18 01:58:40 mstorti Exp $
 
#include <src/debug.h>
#include <cstdio>
#include <time.h>

#include <petsc.h>

Debug debug;

int Debug::stop_f=0;

void Debug::set_signal(int sig) {
  stop_f++;
  if (stop_f>1) {
    signal(SIGINT,orig_handler);
    raise(SIGINT);
  }
}

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
  int stopp;
#define MXTM 100
  char t[MXTM];
  MPI_Comm_rank(comm,&myrank);
  MPI_Allreduce(&stop_f,&stopp,1,MPI_INT,
		MPI_MAX,comm);
  stop_f = stopp;
  if (stop_f) { activate(); stop_f=0;}
  if (myrank==0) {
    if ((active() || active("print"))) {
      tt = time(NULL);
      // t = asctime(localtime(&tt));
      strftime(t,MXTM,"%H:%M:%S",localtime(&tt));
      printf("-- %s -- [%s %10.3f]\n",s,t,chrono.elapsed());
    }
  }
  if (!active()) return;
  int ierr;
  char ans,c;
  ierr = MPI_Barrier(comm);
  assert(ierr==0);
  while (1) {
    if (myrank==0) {
      printf("Command? (cqpd) > ");
      fflush(stdout);
      scanf("%c",&ans);
      printf("\n");
      while (1) { // flush stdin unil newline is found
	scanf("%c",&c);
	if (c=='\n') break;
      }
      // printf("char %c,ord %d\n",ans,ans);
    }
    ierr = MPI_Bcast (&ans, 1, MPI_CHAR, 0,comm);
    assert(ierr==0); 
    if (ans=='q') {
      PetscPrintf(comm,"Quitting by user request...\n");
      PetscFinalize();
      exit(0);
      break;
    } else if (ans=='p') {
      PetscPrintf(comm,"Printing Pid list:\n");
      PetscSynchronizedPrintf(comm,
			      "[%d] Pid: %d\n",myrank,getpid());
      PetscSynchronizedFlush(comm);
    } else if (ans=='d') {
      PetscPrintf(comm,"Deactivate debugging mode...\n");
      deactivate();
      break;
    } else if (ans=='c' || ans=='\n') {
      break;
    } 
  }
  ierr = MPI_Barrier(comm);
  assert(ierr==0);  
}

sighandler_t Debug::orig_handler = NULL;

void Debug::init() {
  orig_handler = signal(SIGINT,&Debug::set_signal);
}
    

Debug::Debug(int active_=0,MPI_Comm comm_=MPI_COMM_WORLD) : 
  comm(comm_) {
  if (active_) {
    activate();
  } else {
    deactivate();
  }
  chrono.start();
}
