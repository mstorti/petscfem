//__INSERT_LICENSE__
//$Id: arglistn.cpp,v 1.2 2001/04/01 01:35:06 mstorti Exp $
 
#include "arglistn.h"

void  OutArg::chunk_pre_operations(int chunk_size,int ndoft) {
  retval = new double[chunk_size*ndoft];
}
