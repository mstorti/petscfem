// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: distmat.h,v 1.9 2001/08/16 03:54:48 mstorti Exp $
#ifndef DISTMAT_H
#define DISTMAT_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This the row of the matrix. 
*/ 
class Row : public  map<int,double> {
  /// define their own members for packing and unpacking
  int size_of_pack() const;
  /// define their own members for packing and unpacking
  void pack(char *&buff) const;
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

void DistMatrix::insert_val(int i,int j,double v) {
  map<int,Row>::iterator I = find(i);
  Row::iterator J;
  if (I == end()) {
    insert(pair<int,Row>(i,Row()));
    I = find(i);
  }
  Row &row = I->second;
  J = row.find(j);
  if (J == row.end()) {
    row.insert(pair<int,double>(j,v));
  } else {
    J->second += v;
  }
}

double DistMatrix::val(int i,int j) {
  Row::iterator J;
  map<int,Row>::iterator I = find(i);
  if (I == end()) return 0.;
  Row &row = I->second;
  J = row.find(j);
  if (J == row.end()) {
    return 0.;
  } else {
    return J->second;
  }
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int DistMat::size_of_pack(map<int,Row>::const_iterator iter) const {
  int n = iter->second.size();
  // size + row number + size*(int+double)
  return (n+2)*sizeof(int)+n*sizeof(double);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistMat::pack(const int &k,const Row &row,char *&buff) const {
  int n = row.size();
  Row::const_iterator iter;
  // number of elements in the row
  //  BufferPack::pack(n,buff);
  BUFFER_PACK(n,buff);
  // the number of the row (key)
  BUFFER_PACK(k,buff);
  // pack the row, each pair (int,double) at a time
  for (iter=row.begin(); iter!=row.end(); iter++) {
    BUFFER_PACK(iter->first,buff);
    BUFFER_PACK(iter->second,buff);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void DistMat::unpack(int &,Row &,const char *&)"
void DistMat::unpack(int &k,Row &row,const char *&buff) {
  int n,j,key;
  double val;
  // Clear the row
  row.clear();
  // unpack the number of elements
  BUFFER_UNPACK(n,buff);
  // unpack the key (number of row)
  BUFFER_UNPACK(k,buff);
  // unpack the pair (key vals) each at a time
  for (j=0; j<n; j++) {
    BUFFER_UNPACK(key,buff);
    BUFFER_UNPACK(val,buff);
    // insert the pair key,val in the row.
    row[key] += val;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
#undef __FUNC__
#define __FUNC__ "int DistMat::processor(const DistMatrix::iterator ) const"
int DistMat::processor(const map<int,Row>::iterator k) const {
  // The number of processor is computed in the `Partitioner' 
  return part->dofpart(k->first);
}
#endif

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
