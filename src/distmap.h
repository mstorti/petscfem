// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: distmap.h,v 1.13 2001/08/10 02:07:14 mstorti Exp $
#ifndef DISTMAP_H
#define DISTMAP_H

#include <map>
#include <vector>
#include <mpi.h>

#include <vecmacros.h>

/** This object class gives information to to which processor belongs
    a dof, a node a processor or whatever. (Currently is only
    implemented for dof's).  
*/
class Partitioner {
public:
  /// Make it pure virtual
  virtual ~Partitioner()=0;
  /* To which processor belongs an element.
     NOT IMPLEMENTED YET
     @param element (input) the element to ask
     @return the processor to which #element# belongs
  */ 
  // virtual int epart(ElementIterator &element) {return 0;};

  /* To which processor belongs a node.
     NOT IMPLEMENTED YET
     @param node (input) the node to ask
     @return the processor to which #node# belongs
  */ 
  virtual int npart(int node) {return 0;};

  /* To which processor belongs a dof.
     NOT IMPLEMENTED YET
     @param dof (input) the dof to ask
     @return the processor to which #dof# belongs
  */ 
  virtual int dofpart(int dof) {return 0;};
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Distributed map class. Elements can be assigned as for a standard
    `map' and, after, a `scatter' operation allows items in the map to
    be passed from one processor to other. The #Partitioner# object
    determines to which processor belongs each dof. 
*/
template <class Key,class Val>
class DistMap : public map<Key,Val> {
 protected:
  /// MPI communicator
  MPI_Comm comm;
  /// This returns the number of processor for a given dof
  Partitioner *part;
  /// size and rank in the comunicator
  int size,myrank;
 public:
  /** Constructor from a communicator
      @param comm_ (input) MPI communicator
      @return a reference to the matrix.
  */ 
  DistMap<Key,Val>(Partitioner *p=NULL,MPI_Comm comm_=MPI_COMM_WORLD);
  /** User defines this function that determine to which processor
      belongs each entry
      @param k (input) iterator to the considered entry. 
      @return the number of processor where these matrix should go. 
  */ 
  int processor(const map<Key,Val>::iterator k) const;
  /** Computes the size of data needed to pack this entry 
      @param k (input) iterator to the entry
      @return the size in bytes of the packed object
   */ 
  int size_of_pack(const map<Key,Val>::iterator k) const;
  /** Packs the entry #(k,v)# in buffer #buff#. This function should
      be defined by the user. 
      @param k (input) key of the entry
      @param v (input) value of the entry
      @param buff (input/output) the position in the buffer where the
      packing is performed
  */ 
  void pack(const Key &k, const Val &v,char *&buff) const;
  /** Does the reverse of #pack#. Given a buffer #buff# recovers the
      corresponding key and val. This function should
      be defined by the user. 
      @param k (output) key of the entry
      @param v (output) value of the entry
      @param buff (input/output) the position in the buffer from where the
      unpacking is performed
  */ 
  void unpack(Key &k,Val &v,const char *& buff);
  /// perform the scatter of elements to its corresponding processor. 
  void scatter();
  /** This function should be defined by the user. Merges a pair key,
      value in the container. 
      @param p (input) the pair to be inserted.
  */ 
  void combine(const pair<Key,Val> &p);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class Key,class Val> DistMap<Key,Val>::
DistMap<Key,Val>(Partitioner *p=NULL,MPI_Comm comm_=MPI_COMM_WORLD) : comm(comm_) {
  // Determine size of the communicator and rank of the processor
  MPI_Comm_size (comm, &size);
  MPI_Comm_rank (comm, &myrank);
  // initialize the partitioner 
  part=p;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// shortcut to access `to_send' as a size x size matrix
#define SEND(p,q) VEC2(to_send,p,q,size)
template <class Key,class Val>
void DistMap<Key,Val>::scatter() {
  map<Key,Val>::iterator iter;
  int *to_send,*to_send_buff,*recv_ok,n_recv_ok,send_ok,
    dest,source;
  pair<Key,Val> p;

  char **send_buff,**send_buff_pos,*recv_buff;
  const char *recv_buff_pos,*recv_buff_pos_end;
  MPI_Request send_rq,recv_rq;
  MPI_Status status;
  int j,k,nsent;
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
  for (iter = begin(); iter != end(); iter++) {
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
  send_buff = new (char *)[size];
  // send_buff_pos:= an array of positions in each of the buffers
  send_buff_pos = new (char *)[size];
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
  for (iter = begin(); iter != end(); iter++) {
    k = processor(iter);
    if (k!=myrank) 
      pack(iter->first,iter->second,send_buff_pos[k]);
  }

  // Erase members that do not belong to this processor.
  for (iter = begin(); iter != end(); iter++) {
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


#if 1 // New scheduling algorithm

  // recv_buff:= buffer for receiving 
  // recv_buff_pos:= positions in the receive buffer
  // recv_buff_pos_end:= end of receive buffer

  // sproc:= mproc:= eproc:= Processes in the lower band (band=0) are
  // s1<= proc< mproc and higher band (band=1) are mproc<= proc <
  // eproc
  int sproc,mproc,eproc,band,stage,s1,s2,nrecv,rank,
    size_here,max_local_size,max_recv_buff_size;

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
  while (1) {
    size_here = eproc-sproc;
    printf("[%d] size here %d\n",myrank,size_here);
    MPI_Allreduce(&size_here,&max_local_size,1,MPI_INT,MPI_MAX,
		  comm);
    if (max_local_size<=1) break;

    if ( size_here> 1) {
      mproc = (sproc+eproc)/2;
      band = (myrank >= mproc);
      printf("[%d] sproc %d, mproc %d, eproc %d\n",
	     myrank,sproc,mproc,eproc);
      for (stage=0; stage<2; stage++) {
	// Range and number of procs in the other band
	if (band==0) {
	  s1=mproc; s2=eproc; 
	} else {
	  s1=sproc; s2=mproc; 
	}
	nrecv=s2-s1;

	if (stage == band) {
	  // Send stage. Send to s1 <= rank < s2
	  for (dest = s1; dest < s2; dest++) {
	    printf("[%d] Sending to %d\n",myrank,dest);
	    MPI_Send(send_buff[dest],SEND(myrank,dest),MPI_CHAR,
		     dest,myrank,comm);
	  }
	} else {
	  // Receive stage
	  printf("[%d] Receive from %d procs.\n",myrank,nrecv);
	  for (rank = 0; rank < nrecv; rank++) {
	    MPI_Recv(recv_buff,max_recv_buff_size,MPI_CHAR,
		     MPI_ANY_SOURCE,MPI_ANY_TAG,
		     comm,&status);
	    source = status.MPI_SOURCE;
	    printf("[%d] received source %d, tag %d\n",
		   myrank,status.MPI_SOURCE,status.MPI_TAG);
	    assert(status.MPI_TAG == source);
	    
	    MPI_Get_count(&status,MPI_CHAR,&nsent);
	    assert(nsent == SEND(source,myrank));

	    recv_buff_pos = recv_buff;
	    recv_buff_pos_end = recv_buff + nsent;
	    while (recv_buff_pos < recv_buff_pos_end ) {
	      unpack(p.first,p.second,recv_buff_pos);
//  	      printf("[%d] unpacking: key %d, val %f\n",
//  		     myrank,p.first,p.second);
	      combine(p);
	    }
	  }
	}
      }

      if (band==0) {
	eproc = mproc;
      } else {
	sproc = mproc;
      }
    }
    // MPI_Barrier(PETSC_COMM_WORLD);
  }

  delete[] recv_buff;

#else // old scheduling algorithm

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
#endif

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
