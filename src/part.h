// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: part.h,v 1.1.2.2 2001/12/24 03:59:56 mstorti Exp $
#ifndef PART_H
#define PART_H

class DofPartitioner {
public:
  virtual int processor(int j) const =0;
  virtual ~DofPartitioner()=0;
};

template <typename ImgValueType>
class Partitioner  {
private:
  typedef pair<int,ImgValueType> ValueType;
public:
  const DofPartitioner *part;
  void processor(const ValueType &p,int &nproc,int *plist) const {
    nproc = 1;
    plist[0] = part->processor(p.first);
  }
  Partitioner<ImgValueType> (const DofPartitioner *pp) 
    : part(pp) { assert(pp!=NULL) ; }
};

#endif
