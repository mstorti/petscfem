//__INSERT_LICENSE__
//$Id: blessf.cpp,v 1.19 2005/01/07 23:39:06 mstorti Exp $

#include <set>

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/pfmat.h>
#include <src/gatherer.h>

#include "advective.h"
#include "nwadvdif.h"
#include "nwadvdifj.h"
#include "burgers.h"
#include "genload.h"
#include "aquifer.h"
#include "stream.h"

//#include "bubbly.h"
#include "bubbly_asm.h"
#include "bubbly_asm_mix.h"
#include "bubbly_full.h"

#include "advec.h"
#include "gasflow.h"
#include "advdfgth.h"
#include "smoke.h"
#include "streamsw1d.h"
#include "nonlres.h"
#include "id.h"
#include "gaschem.h"
#include "./gftwall.h"

#include <time.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "bless_elemset"
void bless_elemset(char *type,Elemset *& elemset) {
  elemset=NULL;
  bless_elemset0(type,elemset);
  if (elemset) return;
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
      //    SET_ELEMSET_TYPE(advdif_swfm1t)
    SET_ELEMSET_TYPE(streamsw1d)
    SET_ELEMSET_TYPE(streamsw1d_abso)

    SET_ELEMSET_TYPE(wall_swfm2t)

    SET_ELEMSET_TYPE(lin_gen_load)

    SET_ELEMSET_TYPE(aquifer)
    SET_ELEMSET_TYPE(stream)
    SET_ELEMSET_TYPE(stream_loss)

    SET_ELEMSET_TYPE(bubbly_full)
    SET_ELEMSET_TYPE(bubbly_bcconv_full)

    SET_ELEMSET_TYPE(bubbly_asm)
    SET_ELEMSET_TYPE(bubbly_bcconv_asm)

    SET_ELEMSET_TYPE(bubbly_asm_mix)
    SET_ELEMSET_TYPE(bubbly_bcconv_asm_mix)

    SET_ELEMSET_TYPE(advec)
    SET_ELEMSET_TYPE(gasflow)
    SET_ELEMSET_TYPE(gasflow_twall)
    SET_ELEMSET_TYPE(gasflow_bcconv)

    SET_ELEMSET_TYPE(flow_rate_integrator)
    SET_ELEMSET_TYPE(id)

    SET_ELEMSET_TYPE(smoke)
    SET_ELEMSET_TYPE(gaschem)
    SET_ELEMSET_TYPE(gaschem_bcconv)
    {
      printf("not known elemset type: \"%s\"\n",type);
      exit(1);
    }
}
