//__INSERT_LICENSE__
/* $Id: bless.cpp,v 1.1 2002/05/04 23:28:15 mstorti Exp $ */

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

