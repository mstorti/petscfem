//__INSERT_LICENSE__
//$Id: dvector.cpp,v 1.1 2003/05/12 02:06:59 mstorti Exp $
 
#include <src/dvector.h>
#include <src/dvector2.h>

// Explicit instantiation for `int' and `double'
template class dvector<double>;
template class dvector<int>;
