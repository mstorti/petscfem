//__INSERT_LICENSE__
/* $Id: blessf.cpp,v 1.19 2003/01/25 17:15:10 mstorti Exp $ */

#include <src/debug.h>
#include <malloc.h>

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/sttfilter.h>
#include <src/pfmat.h>

#include "./nsi_tet.h"
#include "./adaptor.h"
#include "./elast.h"
#include "./qharm.h"
#include "./qharmm.h"
#include <src/gatherer.h>
#include "./nsgath.h"
#include "./embgath.h"
#include "./nssup.h"
#include "./nsikepsrot.h"
#include "./fracstep.h"
#include "./nsid.h"
#include "./mmove.h"
#include "./genload.h"

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
    SET_ELEMSET_TYPE(qharmm)
    SET_ELEMSET_TYPE(ns_id)
    SET_ELEMSET_TYPE(ns_sup)
    SET_ELEMSET_TYPE(ns_sup_g)
    SET_ELEMSET_TYPE(ns_sup_res)
      // SET_ELEMSET_TYPE(elasticity_f)
    SET_ELEMSET_TYPE(elasticity)
    SET_ELEMSET_TYPE(nsi_tet_les_fm2)
    SET_ELEMSET_TYPE(nsi_tet_les_ther)
    SET_ELEMSET_TYPE(nsi_tet_keps)
    SET_ELEMSET_TYPE(nsi_tet_keps_rot)
    SET_ELEMSET_TYPE(fracstep)
    SET_ELEMSET_TYPE(fracstep_fm2)

    SET_ELEMSET_TYPE(bcconv_ns_fm2)
    SET_ELEMSET_TYPE(bcconv_nsther_fm2)
    SET_ELEMSET_TYPE(wall)
    SET_ELEMSET_TYPE(wallke)
    SET_ELEMSET_TYPE(wall_law_res)
    SET_ELEMSET_TYPE(force_integrator)
    SET_ELEMSET_TYPE(visc_force_integrator)
    SET_ELEMSET_TYPE(flow_rate_integrator)
    SET_ELEMSET_TYPE(free_surface_level_integrator)

    SET_ELEMSET_TYPE_ALIAS(mesh_move,mesh_move_eig_anal)
    SET_ELEMSET_TYPE(mesh_move_eig_anal)

    SET_ELEMSET_TYPE(lin_gen_load)
      {
	PETSCFEM_ERROR("not known elemset type: \"%s\"\n",type);
      }
}
