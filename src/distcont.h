// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: distcont.h,v 1.1 2001/08/22 02:25:25 mstorti Exp $
#ifndef DISTCONT_H
#define DISTCONT_H

#include <vector>
//#include <algorithm>

#include <mpi.h>

#include <vecmacros.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Distributed map class. Elements can be assigned as for a standard
    `map' and, after, a `scatter' operation allows items in the map to
    be passed from one processor to other. The #Partitioner# object
    determines to which processor belongs each dof. 
*/
template <typename Container,typename ValueType,typename Partitioner>
class DistCont : public Container {
private:
  int belongs(typename Container::const_iterator k,int *plist) const;
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
  DistCont<Container,
    ValueType,Partitioner>(Partitioner *pp=NULL,
			   MPI_Comm comm_=MPI_COMM_WORLD);
  /** User defines this function that determine to which processor
      belongs each entry
      @param k (input) iterator to the considered entry. 
      @return the number of processor where these matrix should go. 
  */ 
  void processor(const ValueType &p,int &nproc,int *plist) const;
  /** Computes the size of data needed to pack this entry 
      @param k (input) iterator to the entry
      @return the size in bytes of the packed object
  */ 
  int size_of_pack(const ValueType &p) const;
  /** Packs the entry #(k,v)# in buffer #buff#. This function should
      be defined by the user. 
      @param k (input) key of the entry
      @param v (input) value of the entry
      @param buff (input/output) the position in the buffer where the
      packing is performed
  */ 
  void pack(const ValueType &p,char *&buff) const;
  /** Does the reverse of #pack#. Given a buffer #buff# recovers the
      corresponding key and val. This function should
      be defined by the user. 
      @param k (output) key of the entry
      @param v (output) value of the entry
      @param buff (input/output) the position in the buffer from where the
      unpacking is performed
  */ 
  void unpack(ValueType &p,const char *& buff);
  /// perform the scatter of elements to its corresponding processor. 
  void scatter();
  /** This function should be defined by the user. Merges a pair key,
      value in the container. 
      @param p (input) the pair to be inserted.
  */ 
  void combine(const ValueType &p);
#if 0
  class Belongs {
    DistCont *dm;
    int *plist,size,myrank;
  public:
    bool operator() (const ValueType &p) const {
      // bool operator() (typename Container::const_iterator p) const {
      int nproc;
      dm->processor(p,nproc,plist);
      for (int j=0; j<nproc; j++) 
	if (plist[j]==myrank) return false;
      return true;
    }
    void init(DistCont *dm_c,int size_c) {
      dm = dm_c;
      size = size_c;
      plist = new int[size];
    };
    ~Belongs() {delete[] plist;};
  };
#endif
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class Container,typename ValueType,class Partitioner>
DistCont<Container,ValueType,Partitioner>::
DistCont<Container,
  ValueType,
  Partitioner>(Partitioner *pp=NULL, MPI_Comm comm_= MPI_COMM_WORLD) 
    : comm(comm_) {
  // Determine size of the communicator and rank of the processor
  MPI_Comm_size (comm, &size);
  MPI_Comm_rank (comm, &myrank);
  part=pp;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<typename Container,typename ValueType,class Partitioner> 
void DistCont<Container,ValueType,Partitioner>::
processor(const ValueType &p,int &nproc,int *plist) const {
  return part->processor(p,nproc,plist);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<typename Container,typename ValueType,class Partitioner> 
int DistCont<Container,ValueType,Partitioner>::
belongs(typename Container::const_iterator k,int *plist) const {
  int nproc,j;
  processor(*k,nproc,plist);
  for (j=0; j<nproc; j++) 
    if (plist[j] == myrank) return 1;
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// shortcut to access `to_send' as a size x size matrix
#define SEND(p,q) VEC2(to_send,p,q,size)
template <class Container,typename ValueType,class Partitioner>
void DistCont<Container,ValueType,Partitioner>::scatter() {
  typename Container::iterator iter,next;
  int *to_send,*to_send_buff,*recv_ok,n_recv_ok,send_ok,
    dest,source,my_band_start;
  ValueType p;
  char **send_buff,**send_buff_pos,*recv_buff;
  const char *recv_buff_pos,*recv_buff_pos_end;
  MPI_Request send_rq,recv_rq;
  MPI_Status status;
  int j,k,l,nsent;
  int sproc,mproc,eproc,band,stage,s1,s2,nrecv,rank,
    size_here,max_local_size,max_recv_buff_size,jd,nproc,
    eras;
  int *plist = new int[size];

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
    processor(*iter,nproc,plist);
    for (j=0; j<nproc; j++) 
      if (plist[j]!=myrank)
	SEND(myrank,plist[j]) += size_of_pack(*iter);
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
    processor(*iter,nproc,plist);
    for (j=0; j<nproc; j++) {
      k = plist[j];
      if (k!=myrank) pack(*iter,send_buff_pos[k]);
    }
  }

  // Erase members that do not belong to this processor.
  // This implementation may be very slow for non-linked containers. 

  // predicate belongs() tells if the element remains here
#if 1
#if 0
  // Advance until find the first element that remains here
  while (1) {
    iter = begin();
    if (belongs(iter,plist) || iter == end()) break;
    erase(iter);
  }
  if (iter != end()) {
    // This implementation is very careful with respect to not reusing
    // iterators that have been deleted. (That may cause problems). 

    // iter:= keeps on the last element that remains here
    // next:= we compute `next' as `iter++' and check if belongs or
    // not. The we remove it or advence iter. 
    while (1) {
      next = iter;
      next++;
      if (next == end()) break; // Reaches container end
      if (belongs(next,plist)) {
	iter = next;		// advance iterator
      } else {
	erase(next);		// remove item
      }
    }
  }
#endif
#else
  iter = remove_if(begin(),end(),Belongs());
  erase(iter,end());
#endif

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
	    MPI_Send(send_buff[dest],SEND(myrank,dest),MPI_CHAR,
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
	    MPI_Recv(recv_buff,max_recv_buff_size,MPI_CHAR,
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
  delete[] plist;

  // free memory
  delete[] recv_buff;

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
