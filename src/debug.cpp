//__INSERT_LICENSE__
//$Id: debug.cpp,v 1.13 2003/05/04 16:43:50 mstorti Exp $
 
#include <src/debug.h>
#include <sys/resource.h>
#include <sys/vtimes.h>
#include <cstdio>
#include <time.h>

#include <petsc.h>

#include <src/autostr.h>

Debug *GLOBAL_DEBUG=NULL;

int Debug::stop_f=0;

enum PF_MPI_tags { release_barrier };

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

void Debug::release_proc(int proc=0) {
  if (proc<0 || proc>=size) {
    printf("bad procesor number: %d",proc);
    return;
  }
  if (proc==0) {
    for (int k=1; k<size; k++) 
      release_proc(k);
    return;
  }
  if (!flags[proc])
    MPI_Send(&dummy,1,MPI_INT,proc,release_barrier,comm);
  flags[proc]=1;
}

void Debug::trace(const char *s=NULL) {
  time_t tt;
  int stopp;
  char *token;
#define MXTM 100
  char t[MXTM];
  char *wsp = " \n\t"; // whitespace

  if (! was_initialized) {
    MPI_Comm_rank(comm,&myrank);
    MPI_Comm_size(comm,&size);
    was_initialized = 1;
  }

  MPI_Allreduce(&stop_f,&stopp,1,MPI_INT,
		MPI_MAX,comm);
  stop_f = stopp;
  if (stop_f) { activate(); stop_f=0;}
  if (!myrank && (active() || active("print"))) {
    tt = time(NULL);
    // t = asctime(localtime(&tt));
    strftime(t,MXTM,"%H:%M:%S",localtime(&tt));
    printf("-- %s -- [%s %10.3f]\n",s,t,chrono.elapsed());
  }
  if (active("memory_usage")) {
    int ierr;
    char *line;
    int mem,mem_min,mem_max,mem_sum,mem_avrg;
    size_t n=0;
    AutoString file;
    file.sprintf("/proc/%d/status",getpid());
    FILE *fid = fopen(file.str(),"r");
    while(1) {
      assert(getline(&line,&n,fid)!=-1);
      if (sscanf(line,"VmRSS: %d kB",&mem)) break;
    }
    fclose(fid);
    MPI_Allreduce(&mem,&mem_min,1,MPI_INT,MPI_MIN,comm);
    MPI_Allreduce(&mem,&mem_max,1,MPI_INT,MPI_MAX,comm);
    MPI_Allreduce(&mem,&mem_sum,1,MPI_INT,MPI_SUM,comm);
    mem_avrg = int(ceil(double(mem_sum)/double(size)));
    PetscPrintf(PETSC_COMM_WORLD,
		"-- %s -- [Memory usage(kB): min %d, max %d, avrg %d]\n",
		s,mem_min,mem_max,mem_avrg);
  }
  if (!active()) return;
  int ierr,nread,proc;
  char ans,c;
  ierr = MPI_Barrier(comm);
  assert(ierr==0);
  while (1) {
    if (myrank==0) {
      printf("Command? (cqpd) > ");
      fflush(stdout);
      getline (&line,&N,stdin);
      ans = line[0];
      // printf("\n");
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
    } else if (ans=='c') {
      if (myrank==0) {
	flags.clear();
	flags.resize(size,0);
	while (true) {
	  int ntoken=0;
	  token = strtok(&line[1],wsp);
	  while (token!=NULL) {
	    ntoken++;
	    nread = sscanf(token,"%d",&proc);
	    if (nread<=0) break;
	    release_proc(proc);
	    token = strtok(NULL,wsp);
	  }
	  // If no token has been entered then release all processors
	  if (ntoken==0) release_proc();
	  // done:= flags[0] && flags[1] && ... && flags[size-1]
	  int done=1;
	  for (int procc=1; procc<size; procc++) {
	    done = done && flags[procc];
	    if (!done) break;
	  }
	  // If all processors have been released then
	  // release root processor
	  if (done) break;
	  // Expect for more processors to release
	  printf("release procs > ");
	  fflush(stdout);
	  getline (&line,&N,stdin);
	}
      } else {
	MPI_Status stat;
	int flag;
	MPI_Recv(&dummy,1,MPI_INT,0,release_barrier,comm,&stat);
      }
      break;
    } else if (ans=='\n') {
      break;
    } 
  }

  assert(ierr==0);  
}

sighandler_t Debug::orig_handler = NULL;

void Debug::init() {
  orig_handler = signal(SIGINT,&Debug::set_signal);
}

Debug::Debug(int active_=0,MPI_Comm comm_=PETSC_COMM_WORLD) : 
  comm(comm_), was_initialized(0) {
  N = 100;
  line = (char *)malloc(sizeof(char)*N);
  if (active_) {
    activate();
  } else {
    deactivate();
  }
  chrono.start();
}
