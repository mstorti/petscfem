//__INSERT_LICENSE__
//$Id: distcont2.h,v 1.6.10.1 2007/02/19 20:23:56 mstorti Exp $

#ifndef DISTCONT2_H
#define DISTCONT2_H

//#define DEBUG_SPEEDUP

#include <src/distcont.h>
#include <src/utils.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class Container,typename ValueType,class Partitioner>
DistCont<Container,ValueType,Partitioner>::
DistCont(Partitioner *pp, MPI_Comm comm_,iter_mode_t iter_mode_a) 
    : comm(comm_), iter_mode(iter_mode_a) {
  // Determine nprocs of the communicator and rank of the processor
  MPI_Comm_size (comm, &nprocs);
  MPI_Comm_rank (comm, &myrank);
  part=pp;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<typename Container,typename ValueType,class Partitioner> 
int DistCont<Container,ValueType,Partitioner>::
belongs(typename Container::const_iterator k,int *plist) const {
  int nproc,j;
  part->processor(*k,nproc,plist);
  for (j=0; j<nproc; j++) 
    if (plist[j] == myrank) return 1;
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// shortcut to access `to_send' as a nprocs x nprocs matrix
#undef SEND
#define SEND(p,q) VEC2(to_send,p,q,nprocs)
template <class Container,typename ValueType,class Partitioner>
void DistCont<Container,ValueType,Partitioner>::scatter() {
#if DEBUG_SPEEDUP
  HPChrono hpc;
#endif
  typename Container::iterator iter,next;
  int *to_send,*to_send_buff,/**recv_ok,n_recv_ok,send_ok,*/
    dest,source,my_band_start;

  char **send_buff,**send_buff_pos,*recv_buff;
  const char *recv_buff_pos,*recv_buff_pos_end;
  /*MPI_Request send_rq,recv_rq;*/
  MPI_Status status;
  int j,k,/*l,*/nsent;
  int sproc,mproc,eproc,band,stage,s1,s2,nrecv,rank,
    size_here,max_local_size,max_recv_buff_size,jd,nproc;
  int *plist = new int[nprocs];

  // to_send:= `to_send(j,k)' contains the table of how much amount of
  // data has to be sent from processor `j' to processor `k'
  to_send = new int[nprocs*nprocs];
  // auxiliary buffer to perform an `Allreduce'
  to_send_buff = new int[nprocs*nprocs];
  
  // Initialize `to_send' table
  for (j=0; j<nprocs; j++) {
    for (k=0; k<nprocs; k++) {
      SEND(j,k)=0;
    }
  }

  // Compute the table `to_send'
  for (iter = this->begin(); !(iter == this->end()); iter++) {
    part->processor(*iter,nproc,plist);
    for (j=0; j<nproc; j++) 
      if (plist[j]!=myrank)
	SEND(myrank,plist[j]) += size_of_pack(*iter);
  }
  // do not send data to itself
  SEND(myrank,myrank)=0;

  // gather/scatter all `to_send' data 
  MPI_Allreduce(to_send,to_send_buff,nprocs*nprocs,MPI_INT, 
		MPI_SUM,comm);
  // recopy `to_send_buff' to `to_send'
  memcpy(to_send,to_send_buff,nprocs*nprocs*sizeof(int));

#if 0
  // debug: print the `to_send' table
  if (myrank==0) {
    for (j=0; j<nprocs; j++) {
      for (k=0; k<nprocs; k++) {
	if (j != k) printf("%d -> %d: %d\n",j,k,SEND(j,k));
      }
    }
  }
#endif

  // allocate memory for auxiliary vectors
  // send_buff:= An array of buffers for sending
  send_buff = new char *[nprocs];
  // send_buff_pos:= an array of positions in each of the buffers
  send_buff_pos = new char *[nprocs];
  // send_rq:= sendings and receives are non-blocking so that we
  // create a `MPI_Request' object for each of them. 

  // allocate buffers and initialize data
  for (k=0; k<nprocs; k++) {
    // allocate send buffer to proc `k'
    send_buff[k] = new char[SEND(myrank,k)];
    // initialize position 
    send_buff_pos[k] = send_buff[k];
  }

  // Fill send buffers
  for (iter = this->begin(); !(iter == this->end()); iter++) {
    part->processor(*iter,nproc,plist);
    for (j=0; j<nproc; j++) {
      k = plist[j];
      if (k!=myrank) pack(*iter,send_buff_pos[k]);
    }
  }

  if (iter_mode == random_iter_mode) {
    //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    // FOR NON-ASSOCIATIVE CONTAINERS (`erase()' doesn't remove the object,
    // i.e. random-access containers like vectors-deques.)
    for (iter = this->begin(); iter != this->end(); iter++) 
      if (!belongs(iter,plist)) erase(iter);
  } else if (iter_mode == associative_iter_mode) {
    //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
    //FOR ASSOCIATIVE CONTAINERS (`erase()' does remove the object,
    //for instance lists and maps), Advance until find the first
    //element that remains here
    while (1) {
      iter = this->begin();
      if (belongs(iter,plist) || iter == this->end()) break;
      erase(iter);
    }
    if (iter != this->end()) {
      // This implementation is very careful with respect to not reusing
      // iterators that have been deleted. (That may cause problems). 

    // iter:= keeps on the last element that remains here
    // next:= we compute `next' as `iter++' and check if belongs or
    // not. Then we remove it or advance iter. 
      while (1) {
	next = iter;
	next++;
	if (next == this->end()) break; // Reaches container end
	if (belongs(next,plist)) {
	  iter = next;		// advance iterator
	} else {
	  erase(next);		// remove item
	}
      }
    }
  } else assert(0); // bad iter_mode!!

  // Check that all buffers must remain at the end
  for (k=0; k<nprocs; k++) {
    if (k!=myrank && send_buff_pos[k] != send_buff[k]+SEND(myrank,k)) {
      printf("[%d] incorrect position of send buffer to proc %d: \n"
	     "send_buff %p, pos %p, diff %d, to_send %d\n",
	     myrank,k,
	     send_buff[k], send_buff_pos[k],
	     send_buff_pos[k]-send_buff[k], SEND(myrank,k));
    }
  }

  // recv_buff:= buffer for receiving 
  // recv_buff_pos:= positions in the receive buffer
  // recv_buff_pos_end:= end of receive buffer
    
    // Maximu recv buffer nprocs
  max_recv_buff_size=0;
  for (rank=0; rank < nprocs; rank++) {
    if (SEND(rank,myrank) > max_recv_buff_size)
      max_recv_buff_size = SEND(rank,myrank);
  }

  // Alloc recv buffer
  recv_buff = new char[max_recv_buff_size];

  // initially...
  sproc=0;
  eproc=nprocs;

#if DEBUG_SPEEDUP
  hpc.start();
#endif
  // now loop until all groups of processor are of size 1
  while (1) {
    // sproc:= mproc:= eproc:= Processes in the lower band (band=0) are
    // s1<= proc< mproc and higher band (band=1) are mproc<= proc <
    // eproc

    // number of processors in this group
    size_here = eproc-sproc;
    // printf("[%d] size here %d\n",myrank,size_here);
    MPI_Allreduce(&size_here,&max_local_size,1,MPI_INT,MPI_MAX,
		  comm);
    // exit loop if all groups of processor are of size 1
    if (max_local_size<=1) break;

    // Process only if this group has 2 processors at least
    if ( size_here> 1) {
      // Compute middle processor
      mproc = (sproc+eproc)/2;
      // band:= 
      // band = 0 -> I am in the lower band (sproc<= myrank < mproc)
      // band = 1 -> I am in the upper band (mproc<= myrank < eproc)
      band = (myrank >= mproc);
      // printf("[%d] sproc %d, mproc %d, eproc %d\n",
      // myrank,sproc,mproc,eproc);

      // Communication is performed in two stages. In the stage 0
      // the lower band sends to the upper band and the upper band
      // receives. In stage 1 the reverse is performed. 
      for (stage=0; stage<2; stage++) {
	// Range and number of procs in the other band
	if (band==0) {
	  // s1:= s2:= [s1,s2) is the range of processors in the
	  // other band and my_band_start
	  // my_band_start:= where my band starts 
	  s1=mproc; s2=eproc; my_band_start=sproc;
	} else {
	  s1=sproc; s2=mproc; my_band_start=mproc;
	}
	// nrecv:= number of processors in the other band
	nrecv=s2-s1;

	if (stage == band) {
	  // Send stage. Send to s1 <= rank < s2
	  for (jd = 0; jd < nrecv; jd++) {
	    // Shift `dest' to avoid collisions.  In the first
	    // substage my_band_start sends to s1, my_band_start+1
	    // to s1+1 and so on. In the second substage
	    // my_band_start sends to s1+1 and so on...
	    dest = s1 + ((myrank-my_band_start) + jd) % nrecv;
	    // printf("[%d] Sending to %d\n",myrank,dest);
	    int ierr;
	    ierr = MPI_Send(send_buff[dest],SEND(myrank,dest),MPI_CHAR,
			    dest,myrank,comm);
	  }
	} else {
	  // Receive stage
	  // printf("[%d] Receive from %d procs.\n",myrank,nrecv);
	  for (rank = 0; rank < nrecv; rank++) {
	    // We receive packets in any order as there is not
	    // possibility of error since they are tagged with the
	    // source. But we wait for receiving all the packets in
	    // this stage. 
	    int ierr;
	    ierr = MPI_Recv(recv_buff,max_recv_buff_size,MPI_CHAR,
			    MPI_ANY_SOURCE,MPI_ANY_TAG,
			    comm,&status);
	    // Get rank of source 
	    source = status.MPI_SOURCE;
	    // printf("[%d] received source %d, tag %d\n",
	    // myrank,status.MPI_SOURCE,status.MPI_TAG);
	    // assert(status.MPI_TAG == source);

	    nsent = SEND(source,myrank);
#if 1
	    MPI_Get_count(&status,MPI_CHAR,&nsent);
	    assert(nsent == SEND(source,myrank));
#endif
	    // unpack received buffer
	    // recv_buff:= position in buffer
	    recv_buff_pos = recv_buff;
	    // recv_buff_pos_end
	    recv_buff_pos_end = recv_buff + nsent;
	    while (recv_buff_pos < recv_buff_pos_end ) {
	      ValueType p;	// This is to assure that we have a clean
				// `p' each time
	      unpack(p,recv_buff_pos);
	      combine(p);
	    }
	    assert(recv_buff_pos == recv_buff_pos_end);
	  }
	}
      }

      // update the group bounds
      if (band==0) {
	eproc = mproc;
      } else {
	sproc = mproc;
      }
    }
  }
#if DEBUG_SPEEDUP
  PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,
			  "[%d] distcont 1 %f\n",MY_RANK,hpc.elapsed());
  PetscSynchronizedFlush(PETSCFEM_COMM_WORLD);
#endif 

  delete[] plist;

  // free memory
  delete[] recv_buff;

  // Delete all sent and received buffers
  for (k=0; k<nprocs; k++) {
      delete[] send_buff[k];
  }

  // Delete all auxiliary vectors
  delete[] send_buff_pos;
  delete[] send_buff;

  delete[] to_send;
  delete[] to_send_buff;
}

#endif
