// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: distcont.h,v 1.13 2002/08/28 00:48:18 mstorti Exp $
#ifndef DISTCONT_H
#define DISTCONT_H

#include <cstdio>
#include <vector>
#include <mpi.h>
#include <src/vecmacros.h>

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
public:
  /** Iteration is different for non-associative containers (`erase()'
      doesn't remove the object, i.e. random-access containers like
      vectors-deques.) than for associative containers (`erase()' does
      remove the object, for instance lists and maps). (fixme:= this
      should be done better).
  */
  enum iter_mode_t { associative_iter_mode, random_iter_mode };
protected:
  /// MPI communicator
  MPI_Comm comm;
  /// This returns the number of processor for a given dof
  Partitioner *part;
  /// size and rank in the comunicator
  int size,myrank;
  /// Type of iteration mode
  iter_mode_t iter_mode;
public:
  /** Constructor.
      @param part (input) partitioner, defines to which processor belongs each iterator
      @param comm_ (input) MPI communicator
      @param random_iter_mode (input) tells what kind
      of container the class `Container' is. 
  */ 
  DistCont<Container,
    ValueType,Partitioner>(Partitioner *part=NULL,
			   MPI_Comm comm_=MPI_COMM_WORLD,
			   iter_mode_t iter_mode = associative_iter_mode);
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
};

#endif
