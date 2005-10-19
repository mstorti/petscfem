// -*- mode: C++ -*-
//__INSERT_LICENSE__
//$Id: distmap2.h,v 1.12 2005/10/19 17:40:33 mstorti Exp $

#ifndef DISTMAP2_H
#define DISTMAP2_H

#include <src/distmap.h>
//#define DEBUG_SPEEDUP

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class Key,class Val,class Partitioner>
DistMap<Key,Val,Partitioner>::
DistMap<Key,Val,Partitioner>(Partitioner *pp, MPI_Comm comm_) : comm(comm_) {
  // Determine size of the communicator and rank of the processor
  MPI_Comm_size (comm, &size);
  MPI_Comm_rank (comm, &myrank);
  part=pp;
  sched = grouping;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class Key,class Val,class Partitioner> 
int DistMap<Key,Val,Partitioner>::
processor(kv_iterator k) const {
  return part->processor(k);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// shortcut to access `to_send' as a size x size matrix
#define SEND(p,q) VEC2(to_send,p,q,size)
template <class Key,class Val,class Partitioner>
void DistMap<Key,Val,Partitioner>::scatter() {
  HPChrono hpc;
  typename map<Key,Val>::iterator iter;
  // kv_iterator iter;
  int *to_send,*to_send_buff,*recv_ok,n_recv_ok,send_ok,
    dest,source,my_band_start;
  pair<Key,Val> p;

  char **send_buff,**send_buff_pos,*recv_buff;
  const char *recv_buff_pos,*recv_buff_pos_end;
  MPI_Request send_rq,recv_rq;
  MPI_Status status;
  int j,k,nsent;
  int sproc,mproc,eproc,band,stage,s1,s2,nrecv,rank,
    size_here,max_local_size,max_recv_buff_size,jd;

  // to_send:= `to_send(j,k)' contains the table of how much amount of
  // data has to be sent from processor `j' to processor `k'
  to_send = new int[size*size];
  // auxiliary buffer to perform an `Allreduce'
  to_send_buff = new int[size*size];
  
  // Initialize `to_send' table
  for (j=0; j<size; j++) {
    for (k=0; k<size; k++) {
      SEND(j,k)=0;
    }
  }

  // Compute the table `to_send'
  for (iter = this->begin(); iter != this->end(); iter++) {
    k = processor(iter);
    assert(k>=0);
    assert(k<size);
    SEND(myrank,k) += size_of_pack(iter);
  }
  // do not send data to itself
  SEND(myrank,myrank)=0;

  // gather/scatter all `to_send' data 
  MPI_Allreduce(to_send,to_send_buff,size*size,MPI_INT, 
		MPI_SUM,comm);
  // recopy `to_send_buff' to `to_send'
  memcpy(to_send,to_send_buff,size*size*sizeof(int));

#if 0
  // debug: print the `to_send' table
  if (myrank==0) {
    for (j=0; j<size; j++) {
      for (k=0; k<size; k++) {
	if (j != k) printf("%d -> %d: %d\n",j,k,SEND(j,k));
      }
    }
  }
#endif

  // allocate memory for auxiliary vectors
  // send_buff:= An array of buffers for sending
  send_buff = new char *[size];
  // send_buff_pos:= an array of positions in each of the buffers
  send_buff_pos = new char *[size];
  // send_rq:= sendings and receives are non-blocking so that we
  // create a `MPI_Request' object for each of them. 

  // allocate buffers and initialize data
  for (k=0; k<size; k++) {
    // allocate send buffer to proc `k'
    send_buff[k] = new char[SEND(myrank,k)];
    // initialize position 
    send_buff_pos[k] = send_buff[k];
  }

  // Fill send buffers
  for (iter = this->begin(); iter != this->end(); iter++) {
    k = processor(iter);
    if (k!=myrank) 
      pack(iter->first,iter->second,send_buff_pos[k]);
  }

  // Erase members that do not belong to this processor.
  for (iter = this->begin(); iter != this->end(); iter++) {
    k = processor(iter);
    if (k!=myrank) erase(iter);
  }

  // Check that all buffers must remain at the end
  for (k=0; k<size; k++) {
    if (k!=myrank && send_buff_pos[k] != send_buff[k]+SEND(myrank,k)) {
      printf("[%d] incorrect position of send buffer to proc %d: \n"
	     "send_buff %p, pos %p, diff %d, to_send %d\n",
	     myrank,k,
	     send_buff[k], send_buff_pos[k],
	     send_buff_pos[k]-send_buff[k], SEND(myrank,k));
    }
  }


  if (sched == grouping) { // New scheduling algorithm

    // recv_buff:= buffer for receiving 
    // recv_buff_pos:= positions in the receive buffer
    // recv_buff_pos_end:= end of receive buffer
    
    // Maximu recv buffer size
    max_recv_buff_size=0;
    for (rank=0; rank < size; rank++) {
      if (SEND(rank,myrank) > max_recv_buff_size)
	max_recv_buff_size = SEND(rank,myrank);
    }

    // Alloc recv buffer
    recv_buff = new char[max_recv_buff_size];

    // initially...
    sproc=0;
    eproc=size;
    // now loop until all groups of processor are of size 1
    HPChrono hpc2;
    double tsend,trecv;
    hpc2.start();
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
	      hpc.start();
	      MPI_Send(send_buff[dest],SEND(myrank,dest),MPI_CHAR,
		       dest,myrank,comm);
	      tsend += hpc.elapsed();
	    }
	  } else {
	    // Receive stage
	    // printf("[%d] Receive from %d procs.\n",myrank,nrecv);
	    for (rank = 0; rank < nrecv; rank++) {
	      // We receive packets in any order as there is not
	      // possibility of error since they are tagged with the
	      // source. But we wait for receiving all the packets in
	      // this stage. 
	      hpc.start();
	      MPI_Recv(recv_buff,max_recv_buff_size,MPI_CHAR,
		       MPI_ANY_SOURCE,MPI_ANY_TAG,
		       comm,&status);
	      trecv += hpc.elapsed();
	      // Get rank of source 
	      source = status.MPI_SOURCE;
	      // printf("[%d] received source %d, tag %d\n",
	      // myrank,status.MPI_SOURCE,status.MPI_TAG);
	      // assert(status.MPI_TAG == source);

	      nsent = SEND(source,myrank);
#if 0	    
	      MPI_Get_count(&status,MPI_CHAR,&nsent);
	      assert(nsent == SEND(source,myrank));
#endif
	      // unpack received buffer
	      // recv_buff:= position in buffer
	      recv_buff_pos = recv_buff;
	      // recv_buff_pos_end
	      recv_buff_pos_end = recv_buff + nsent;
	      while (recv_buff_pos < recv_buff_pos_end ) {
		unpack(p.first,p.second,recv_buff_pos);
		//  	      printf("[%d] unpacking: key %d, val %f\n",
		//  		     myrank,p.first,p.second);
		combine(p);
	      }
	      // assert(recv_buff_pos == recv_buff_pos_end);
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
#ifdef DEBUG_SPEEDUP
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			    "[%d] in distmap.h, send %f, recv %f, total %f\n",
			    MY_RANK,tsend,trecv,hpc2.elapsed());
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
#endif 
    // free memory
    delete[] recv_buff;

  } else if (sched == rotate_all) {

    // recv_buff:= buffer for receiving 
    // recv_buff_pos:= positions in the receive buffer
    // recv_buff_pos_end:= end of receive buffer
    for (k=1; k<size; k++) { 
      dest = (myrank+k) % size;
      source = (myrank-k+size) % size;

      // allocate recv buffer to from proc  `k'
      recv_buff = new char[SEND(source,myrank)];
      // initialize position 
      recv_buff_pos = recv_buff;
      if (myrank!=0) {
	MPI_Send(send_buff[dest],SEND(myrank,dest),MPI_CHAR,
		 dest,myrank,comm);
	MPI_Recv(recv_buff,SEND(source,myrank),MPI_CHAR,source,source,
		 comm,&status);
      } else {
	MPI_Recv(recv_buff,SEND(source,myrank),MPI_CHAR,source,source,
		 comm,&status);
	MPI_Send(send_buff[dest],SEND(myrank,dest),MPI_CHAR,
		 dest,myrank,comm);
      }    

      MPI_Get_count(&status,MPI_CHAR,&nsent);
      // printf("[%d] %d received from %d\n",myrank,nsent,source);
      if (nsent!=SEND(source,myrank)) 
	printf("[%d] Didn't receive expected amount of data\n"
	       "expected %d, received  %d\n",
	       myrank,SEND(k,myrank),nsent);
      recv_buff_pos_end = recv_buff + SEND(source,myrank);
      while (recv_buff_pos < recv_buff_pos_end ) {
	unpack(p.first,p.second,recv_buff_pos);
	//        PetscPrintf(PETSC_COMM_WORLD,"unpacking: key %d, val %f\n",
	//  		  p.first,p.second);
	combine(p);
      }
      delete[] recv_buff;
    }
  } else {
    assert(0);
  }

  // Delete all sent and received buffers
  for (k=0; k<size; k++) {
      delete[] send_buff[k];
  }

  // Delete all auxiliary vectors
  delete[] send_buff_pos;
  delete[] send_buff;

  delete[] to_send;
  delete[] to_send_buff;
}

#endif
