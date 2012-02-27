// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: qharmm.h,v 1.6 2006/10/23 02:43:18 mstorti Exp $
#ifndef PETSCFEM_LUBRICATION_H
#define PETSCFEM_LUBRICATION_H

#include <src/fm2temp.h>
#include <src/gatherer.h>
#include <src/fastmat2.h>
#include <src/secant.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Adaptation of quasi-harmonic element for lubrication theory
/// (Reynolds equation)
class lubrication : public adaptor_pg {
private:
  double Dt,rec_Dt;
  FastMat2 cond,C,Cp,x_ref,G;
  FastMat2Tmp tmp;
  double 
  L,R,                // Length and radious of bearing
  rho,viscosity,nu,c, // density,viscosity,kin. viscosity, clearance
    Omega0,Omega1, // rotation velocities for stator and rotor
    e0x,e0y,e1x,e1y; // displacements of center for stator and rotor
public:
  void lub_init();
  void elemset_init();
  void pg_connector(const FastMat2 &xpg,
		    const FastMat2 &state_old_pg,
		    const FastMat2 &grad_state_old_pg,
		    const FastMat2 &state_new_pg,
		    const FastMat2 &grad_state_new_pg,
		    FastMat2 &res_pg,FastMat2 &mat_pg);
  void set_pg_values(vector<double> &pg_values,FastMat2 &u,
                     FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
                     double wpgdet,double time);
};

extern lubrication *lubrication_p;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Computes the force on the rotor
class lub_force_integrator : public gatherer {
private:
public:
  /// perform several checks and initialization
  void init();
  /// add volume
  void set_pg_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
		     double wpgdet,double time);
};

#endif
