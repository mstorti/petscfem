//__INSERT_LICENSE__
//$Id: blessf.cpp,v 1.3 2002/02/15 12:47:37 mstorti Exp $

#include <set>

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/pfmat.h>

#include "advective.h"
#include "nwadvdif.h"
#include "nwadvdifj.h"
#include "burgers.h"
#include "genload.h"
#include "aquifer.h"
#include "stream.h"
#include "bubbly.h"

#include <time.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "bless_elemset"
void bless_elemset(char *type,Elemset *& elemset) {
    // General linear advective-diffusive system
    SET_ELEMSET_TYPE(advdif_advecfm2)
    SET_ELEMSET_TYPE(bcconv_adv_advecfm2)
    // new version
    SET_ELEMSET_TYPE(newadvdif_advecfm2)
    SET_ELEMSET_TYPE(newbcconv_advecfm2)
    // Burger's eq.
    SET_ELEMSET_TYPE(bcconv_adv_burgers)
    SET_ELEMSET_TYPE(advdif_burgers)
    // new version
    SET_ELEMSET_TYPE(newadvdif_burgers)
    SET_ELEMSET_TYPE(newbcconv_burgers)
    // Turbulent shallow water
    SET_ELEMSET_TYPE(bcconv_adv_swfm2t)
    SET_ELEMSET_TYPE(advdif_swfm2t)
    SET_ELEMSET_TYPE(wall_swfm2t)

    SET_ELEMSET_TYPE(lin_gen_load)

    SET_ELEMSET_TYPE(aquifer)
    SET_ELEMSET_TYPE(stream)
    SET_ELEMSET_TYPE(stream_loss)

    SET_ELEMSET_TYPE(bubbly)
    {
      printf("not known elemset \"type\": %s\n",type);
      exit(1);
    }
}

