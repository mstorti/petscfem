//__INSERT_LICENSE__
/* $Id: blessf.cpp,v 1.1 2002/05/04 23:31:18 mstorti Exp $ */

#include <src/debug.h>
#include <malloc.h>

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/sttfilter.h>
#include <src/pfmat.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "elast.h"
#include "qharm.h"
#include "gatherer.h"
#include "nssup.h"
#include "nsikepsrot.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "bless_elemset"
void bless_elemset(char *type,Elemset *& elemset) {
  //  SET_ELEMSET_TYPE(internal)
  //    SET_ELEMSET_TYPE(fracstep) Por ahora lo desactivamos hasta que
  // hagamos la interfase
  // SET_ELEMSET_TYPE(nsi_tet)
  //  SET_ELEMSET_TYPE(nsi_tet_les)
    SET_ELEMSET_TYPE(qharm)
    SET_ELEMSET_TYPE(ns_id)
    SET_ELEMSET_TYPE(ns_sup)
    SET_ELEMSET_TYPE(ns_sup_res)
      // SET_ELEMSET_TYPE(elasticity_f)
    SET_ELEMSET_TYPE(elasticity)
    SET_ELEMSET_TYPE(nsi_tet_les_fm2)
    SET_ELEMSET_TYPE(nsi_tet_les_ther)
    SET_ELEMSET_TYPE(nsi_tet_keps)
    SET_ELEMSET_TYPE(nsi_tet_keps_rot)
    SET_ELEMSET_TYPE(nsi_rot)
    SET_ELEMSET_TYPE(bcconv_ns_fm2)
    SET_ELEMSET_TYPE(bcconv_nsther_fm2)
    SET_ELEMSET_TYPE(wall)
    SET_ELEMSET_TYPE(wallke)
    SET_ELEMSET_TYPE(wall_law_res)
    SET_ELEMSET_TYPE(force_integrator)
    SET_ELEMSET_TYPE(flow_rate_integrator)
	{
	printf("not known elemset type: \"%s\"\n",type);
	exit(1);
	}
}

