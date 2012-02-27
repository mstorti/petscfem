//__INSERT_LICENSE__
/* $Id$ */

#include <src/debug.h>
#include <malloc.h>

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/sttfilter.h>
#include <src/pfmat.h>
#include <src/pfobject.h>
#include <src/gatherer.h>
#include <src/srfgath.h>

#include "./nsi_tet.h"

#include "./adaptor.h"
#include "./elast.h"
#include "./elast2.h"
#include "./elastld.h"
#include "./qharm.h"
#include "./errestim.h"
#include "./qharmm.h"
#include "./lubrication.h"
#include "./nsgath.h"
#include "./embgath.h"
#include "./nssup.h"
#include "./nsikepsrot.h"
#include "./fracstep.h"
#include "./nsid.h"
#include "./mmoveopt.h"
#include "./mmoveopt2.h"
#include "./mmoveopt3.h"
#include "./mmove.h"
#include "./mmove2.h"
#include "./genload.h"
#include "./flowrev.h"
#include "./invcoupl.h"
#include "./nullvort.h"
#include "./interplns.h"
#include "./condwall.h"
#include "./condwallpen.h"
#include "./bubblyqint.h"
#include "./bubblyqint.h"
#include "./nsitetlesf.h"
#include "./nsitetlesls.h"
#include "./nsitetlesfbf.h"
#include "./nsitetlesfd.h"
#include "./truss.h"
#include "./nodeload.h"
#include "./poiboltz.h"
#include "./electrophoresisM2.h"
#include "./renorm.h"
#include "./renorm2.h"
#include "./renorm3.h"
#include "./poisson.h"
#include "./charge_cons.h"
#include "./pot_grad.h"
#include "./electrophoresisM.h"
#include "./electrophoresis_mov.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "bless_elemset_ns"
void bless_elemset_ns(char *type,Elemset *& elemset) {
  //  SET_ELEMSET_TYPE(internal)
  //    SET_ELEMSET_TYPE(fracstep) Por ahora lo desactivamos hasta que
  // hagamos la interfase
  // SET_ELEMSET_TYPE(nsi_tet)
  //  SET_ELEMSET_TYPE(nsi_tet_les)
    SET_ELEMSET_TYPE(qharm)
    SET_ELEMSET_TYPE(error_estimator)
    SET_ELEMSET_TYPE(qharmm)

    SET_ELEMSET_TYPE(lubrication)
    SET_ELEMSET_TYPE(lub_force_integrator)

    SET_ELEMSET_TYPE(ns_id)
    SET_ELEMSET_TYPE(nodeload)
    SET_ELEMSET_TYPE(ns_sup)
    SET_ELEMSET_TYPE(ns_sup_g)
    SET_ELEMSET_TYPE(ns_sup_res)
      // SET_ELEMSET_TYPE(elasticity_f)
    SET_ELEMSET_TYPE(elasticity)
    SET_ELEMSET_TYPE(elasticity2)
    SET_ELEMSET_TYPE(ld_elasticity)
    SET_ELEMSET_TYPE(ld_elasticity_load)
    SET_ELEMSET_TYPE_ALIAS(ld_elasticity_df,ld_elasticity)

    SET_ELEMSET_TYPE(nsi_tet_les_fm2)
    SET_ELEMSET_TYPE(nsi_tet_les_ls)
    SET_ELEMSET_TYPE(nsi_tet_les_full)
    SET_ELEMSET_TYPE(nsi_tet_les_full_ctx2)
    SET_ELEMSET_TYPE(nsi_tet_les_full_darcy)
    SET_ELEMSET_TYPE(nsi_tet_les_full_bf)
    SET_ELEMSET_TYPE(nsi_tet_les_comp)
    SET_ELEMSET_TYPE(nsi_tet_les_ther)
    SET_ELEMSET_TYPE(nsi_tet_les_asm)
    SET_ELEMSET_TYPE(bubbly_flow_rate_integrator)
    SET_ELEMSET_TYPE(nsi_tet_keps)
    SET_ELEMSET_TYPE(nsi_tet_keps_rot)
    SET_ELEMSET_TYPE(nsi_tet_asm)
    SET_ELEMSET_TYPE(nsi_tet_asm_avgvol)
    SET_ELEMSET_TYPE(ns_gasflow)

    SET_ELEMSET_TYPE(fracstep)
    SET_ELEMSET_TYPE(fracstep_fm2)
    SET_ELEMSET_TYPE(bcconv_fstep_fm2)
    SET_ELEMSET_TYPE(fracstep_fm2_cw)

    SET_ELEMSET_TYPE(bcconv_ns_fm2)
    SET_ELEMSET_TYPE(bcconv_nsi_tet_asm)
    SET_ELEMSET_TYPE(bcconv_nsi_tet_asm_avgvol)
    SET_ELEMSET_TYPE(bcconv_nsther_fm2)
    SET_ELEMSET_TYPE(bcconv_nsasm_fm2)
    SET_ELEMSET_TYPE(bcconv_ns_gasflow)

    SET_ELEMSET_TYPE(wall)
    SET_ELEMSET_TYPE(wallke)
    SET_ELEMSET_TYPE(wall_law_res)
    SET_ELEMSET_TYPE(force_integrator)
    SET_ELEMSET_TYPE(visc_force_integrator)
    SET_ELEMSET_TYPE(flow_rate_integrator)
    SET_ELEMSET_TYPE(free_surface_level_integrator)
    SET_ELEMSET_TYPE(volume_integrator)
    SET_ELEMSET_TYPE(elast_energy_integrator)

    SET_ELEMSET_TYPE_ALIAS(mesh_move,mesh_move_eig_anal)
    SET_ELEMSET_TYPE(mesh_move_eig_anal)
    SET_ELEMSET_TYPE(mesh_move2)
    SET_ELEMSET_TYPE(mesh_move_opt)
    SET_ELEMSET_TYPE(mesh_move_opt2)
    SET_ELEMSET_TYPE(mesh_move_opt3)
    SET_ELEMSET_TYPE(truss)

    SET_ELEMSET_TYPE(renorm)
    SET_ELEMSET_TYPE(renorm2)
    SET_ELEMSET_TYPE(renorm3)

    SET_ELEMSET_TYPE(inviscid_coupling)

    SET_ELEMSET_TYPE(lin_gen_load)
    SET_ELEMSET_TYPE(flow_reversal)

    SET_ELEMSET_TYPE(cond_wall)
    SET_ELEMSET_TYPE(cond_wall_pen)
    SET_ELEMSET_TYPE(cond_wall_lm)

    SET_ELEMSET_TYPE(dl_penalize)

    SET_ELEMSET_TYPE_ALIAS(interpolation,interpolation_ns)

    SET_ELEMSET_TYPE(poisson_boltzmann)

    SET_ELEMSET_TYPE(poisson)
    
    SET_ELEMSET_TYPE(electrophoresisM2)

    SET_ELEMSET_TYPE(electrophoresisM)

    SET_ELEMSET_TYPE(electrophoresis_mov)

    SET_ELEMSET_TYPE(charge_cons)

    SET_ELEMSET_TYPE(pot_grad)

    {
      elemset=NULL;
    }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
BasicObject *BasicObject_ns_factory(string &type) {
  if (0) {} // tricky!!
  else if (type=="null_vort") return new null_vort_bo;
  else return NULL;
}
