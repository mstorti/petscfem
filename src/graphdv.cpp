//__INSERT_LICENSE__
//$Id: graphdv.cpp,v 1.1 2002/07/22 12:08:32 mstorti Exp $

#include <src/graphdv.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void graphdv::resync() { 
  if (modif) {
    // printf("resyncing at size: %d\n",da.size());
    da.sort();
    da.shrink(da.remove_unique());
    modif = 0;
    // printf("end resync.\n",da.size());
  }
}

