// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: distmap.h,v 1.2 2001/07/30 03:42:27 mstorti Exp $
#ifndef DISTMAP_H
#define DISTMAP_H

#include <map>
#include <vector>
#include <mpi.h>

#include <vecmacros.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Distributed map class. Elements can be assigned as for a standard
    `map' and, after, a `scatter' operation allows items in the map to
    be passed from one processor to other. 
*/
template <class Key,class Val>
class DistMap : public map<Key,Val> {
 private:
  /// MPI communicator
  MPI_Comm comm;
  /// size and rank in the comunicator
  int size,myrank;
 public:
  /** Constructor from a communicator
      @param comm_ (input) MPI communicator
      @return a reference to the matrix.
  */ 
  DistMap<Key,Val>(MPI_Comm comm_=MPI_COMM_WORLD);
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
  void pack(const Key &k, const Val &v,char **buff) const;
  /** Does the reverse of #pack#. Given a buffer #buff# recovers the
      corresponding key and val. This function should
      be defined by the user. 
      @param k (output) key of the entry
      @param v (output) value of the entry
      @param buff (input/output) the position in the buffer from where the
      unpacking is performed
  */ 
  void unpack(Key &k,Val &v,const char *buff);
  /// perform the scatter of elements to its corresponding processor. 
  void scatter();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class Key,class Val> DistMap<Key,Val>::
DistMap<Key,Val>(MPI_Comm comm_=MPI_COMM_WORLD) : comm(comm_) {
  MPI_Comm_size (comm, &size);
  MPI_Comm_rank (comm, &myrank);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// shortcut to access `to_send' as a size x size matrix
#define SEND(p,q) VEC2(to_send,p,q,size)
template <class Key,class Val>
void DistMap<Key,Val>::scatter() {
  map<Key,Val>::iterator iter;
  int *to_send,*to_send_buff,*recv_ok,n_recv_ok,send_ok;
  char **send_buff,**send_buff_pos,**recv_buff,**recv_buff_pos;
  MPI_Request *send_rq,*recv_rq;
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

#if 1
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
  send_rq = new MPI_Request[size];

  // recv_buff:= An array of buffers for receiving
  recv_buff = new (char *)[size];
  // recv_buff_pos:= an array of positions in each of the buffers
  recv_buff_pos = new (char *)[size];
  // Request objects for receiving
  recv_rq = new MPI_Request[size];
  // recv_ok:= flags whether the receive from processor `k' to this
  // has been performed
  recv_ok = new int[size];

  // allocate buffers and initialize data
  for (k=0; k<size; k++) {
    // allocate send buffer to proc `k'
    send_buff[k] = new char[SEND(myrank,k)];
    // initialize position 
    send_buff_pos[k] = send_buff[k];
    // allocate recv buffer to from proc  `k'
    recv_buff[k] = new char[SEND(k,myrank)];
    // initialize position 
    recv_buff_pos[k] = recv_buff[k];
    // initialize flags
    recv_ok[k]=0;
  }

  // Fill send buffers
  for (iter = begin(); iter != end(); iter++) {
    k = processor(iter);
    if (k!=myrank) 
      pack(iter->first,iter->second,&send_buff_pos[k]);
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

  // MPI_Send(&a(row_start,0),nrows_sent*N,MPI_DOUBLE,
  // proc,row_start,MPI_COMM_WORLD);

  // Send all buffers non-blocking
  for (k=0; k<size; k++) {
    if (k!=myrank) {
      printf("[%d] sending %d\n",myrank,k);
      MPI_Isend(send_buff[k],SEND(myrank,k),MPI_CHAR,
	       k,myrank,comm,&send_rq[k]);
    }
  }

  // receive all buffers non-blocking
  for (k=0; k<size; k++) {
    if (k!=myrank) {
      printf("[%d] receiving %d\n",myrank,k);
      MPI_Irecv(recv_buff[k],SEND(k,myrank),MPI_CHAR,k,myrank,
		comm,&recv_rq[k]);
    }
  }

  // Loop until all buffers have been received
  k=0;
  // n_recv_ok:= number of buffers that have been received so far 
  n_recv_ok = 0;
  while (1) {
    // test only non-received buffers
    if (k==myrank || recv_ok[k]) continue;
    printf("[%d] testing %d\n",myrank,k);
    // test unreceived buffer
    MPI_Test(&recv_rq[k],&recv_ok[k],&status);
    if (recv_ok[k]) {
      // Message has arrived
      n_recv_ok++;
      // check that the length of the received message is the expected 
      MPI_Get_count(&status,MPI_CHAR,&nsent);
      printf("[%d] %d received from %d\n",myrank,nsent,k);
      if (nsent!=SEND(myrank,k)) 
	printf("Didn't receive expected amount of data\n");
      // Exit if all buffers have been received
      if (n_recv_ok == size - 1) break;
    }
    k++;
    // k goes cyclically until all messages have been received
    if (k==size) k = 0;
  }

  // Check that all data buffers sent have been sent
  for (k=0; k<size; k++) {
    if (k!=myrank) {
      MPI_Test(&send_rq[k],&send_ok,&status);
      printf("[%d] sent to %d OK ? %d\n",myrank,k,send_ok);
      if (!send_ok) 
	printf("[%d] not sent OK data to %d\n",myrank,k);
    }
  }

  // Delete all sent and received buffers
  for (k=0; k<size; k++) {
      delete[] send_buff[k];
      delete[] recv_buff[k];
  }

  // Delete all auxiliary vectors
  delete[] recv_ok;

  delete[] send_buff_pos;
  delete[] send_buff;
  delete[] send_rq;

  delete[] recv_buff_pos;
  delete[] recv_buff;
  delete[] recv_rq;

  delete[] to_send;
  delete[] to_send_buff;
}

#endif
