// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: part.h,v 1.1.2.1 2001/12/17 00:03:57 mstorti Exp $
#ifndef PART_H
#define PART_H

class DofPartitioner {
public:
  virtual int processor(int j)=0;
};

#endif
