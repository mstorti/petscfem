// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: distmat.h,v 1.10 2001/08/16 18:24:46 mstorti Exp $
#ifndef DISTMAT_H
#define DISTMAT_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This the row of the matrix. 
*/ 
class Row : public  map<int,double> {
private:
  /// define their own members for packing and unpacking
  int size_of_pack() const;
  /// define their own members for packing and unpacking
  void pack(char *&buff) const;
public:
  /// print the row
  void print() const;
};

class IntRowPartitioner {
public:
  virtual ~IntRowPartitioner()=0;
  virtual int processor(map<int,Row>::iterator k)=0;
};

typedef DistMap<int,Row,IntRowPartitioner> DistMat;

/** A distributed map<int,int,double> class
 */ 
class DistMatrix : public DistMat {
public:
  /// Specific ``insert'' routine. 
  void insert_val(int i,int j,double v);
  /// Specific function for retrieving values. 
  double val(int i,int j);
  DistMatrix(IntRowPartitioner *pp=NULL,MPI_Comm comm_=MPI_COMM_WORLD) :
    DistMat(pp,comm_) {};
};


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
#undef __FUNC__
#define __FUNC__ "int DistMat::processor(const DistMatrix::iterator ) const"
int DistMat::processor(const map<int,Row>::iterator k) const {
  // The number of processor is computed in the `Partitioner' 
  return part->dofpart(k->first);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistMat::combine(const pair<int,Row> &p) {
  // Insert the pair (int,row) in the matrix
  int n,j;
  map<int,Row>::iterator iter = find(p.first);
  Row::iterator r;
  Row::const_iterator q;
  if (iter == end()) {
    // insert(p); // this doesn't work
    // If entry doesn't exist for this key, then insert a new row for
    // this key 
    (*this)[p.first] = p.second;
  } else {
    // merge the elements of this row to the existing ones
    // existing row
    Row &oldr = iter->second;
    // new row
    const Row &newr = p.second;
    // new element in the row
    for (q = newr.begin(); q!=newr.end(); q++) {
      // r iterator to that key in the old row
      r = oldr.find(q->first);
      if (r == oldr.end()) {
	// the col is not in the row
	// insert a new entry
	oldr.insert(*q);
      } else {
	// add to the existing 
	r->second += q->second;
      }
    }
  }
}
#endif

#endif
