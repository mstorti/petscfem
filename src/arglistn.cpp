//__INSERT_LICENSE__
//$Id: arglistn.cpp,v 1.4 2001/05/30 18:21:53 mstorti Exp $
 
#include "arglistn.h"

void  OutArg::chunk_pre_operations(int chunk_size,int ndoft) {
  retval = new double[chunk_size*ndoft];
}
