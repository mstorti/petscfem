/*__INSERT_LICENSE__*/
//$Id: sttfilter.cpp,v 1.7 2001/05/30 03:58:50 mstorti Exp $

#include "sttfilter.h"


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void LowPass::update(Time)" 
void LowPass::update(Time time) {
  if (input->step() <= step()) input->update(time); 
  i_state.scale(gamma).axpy(1-gamma,input->state());
  i_state.set_time(time);
  Filter::update(time);
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Mixer & Mixer::add_input(Filter &input,double g)" 
Mixer & Mixer::add_input(Filter &input,double g) {
  filter_l.push_back(&input);
  gain_l.push_back(g);
  return *this;
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void Mixer::update(Time time)" 
void Mixer::update(Time time) {
  i_state.set_cnst(0.);
  for (int j=0; j<filter_l.size(); j++) {
    Filter *f = filter_l[j];
    if (f->step() <= step()) f->update(time); 
    // f->update(time);
    i_state.axpy(gain_l[j],f->state());
    i_state.set_time(time);
  }
  Filter::update(time);
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "LPFilterGroup::LPFilterGroup(TextHashTable *,State &,double)"
LPFilterGroup::LPFilterGroup(TextHashTable *thash,
			     State &x,double Dt) : input(x), nalpha(0) {

  //o _T: double[2*nalpha]
  //  _N: low_pass_filter _D: no filter (nalpha=0)  _DOC: 
  //  Enter pair of values #gamma_1 n_1 gamma_2 n_2# ...
  //  so that the relaxation time are $\tau_j = 1./gamma_j$ and $n_j$
  //  is the corresponding order. 
  //  FIXME:= TO BE DOCUMENTED LATER
  //  _END
  const char *entry;
  vector<double> array;
  VOID_IT(array);
  thash->get_entry("low_pass_filter",entry);
  if (entry) read_double_array(array,entry);
  nalpha = array.size();
  if (int(nalpha/2)*2 != nalpha) {
    PetscPrintf(PETSC_COMM_WORLD,
		"Number of entries in \"low_pass_filter\" should be even\n"
		"Entered %d values\n",nalpha);
    exit(1);
  }
  nalpha /= 2;
  gamma_v.resize(nalpha);
  n_v.resize(nalpha);
  for (int j=0; j<nalpha; j++) {
    double aa = array[2*j];
    double a = array[2*j+1];
    printf("aa %f, a %f\n",aa,a);
    if (aa <= 0.) {
      PetscPrintf(PETSC_COMM_WORLD,
		  "Odd entries in \"low_pass_filter\" should be "
		  "positive reals. Entered %f\n",aa);
      exit(1);
    }
    gamma_v[j] = exp(-aa*Dt);
    int n = int(a);
    if (a != double(n) || n<=0) {
      PetscPrintf(PETSC_COMM_WORLD,
		  "Even entries in \"low_pass_filter\" should be "
		  "positive integers\n"
		  "gamma_1 n_1 gamma_2 n_2 ... Entered %f\n",
		  "Entered %d values\n",a);
      exit(1);
    }
    n_v[j] = n;
  }
  
  lp_filters.resize(nalpha);
  // Create filters1
  for (int j=0; j<nalpha; j++) {
    lp_filters[j].resize(n_v[j]);
    lp_filters[j][0] = new LowPass (gamma_v[j],input,x);
    for (int k=1; k<n_v[j]; k++) {
      lp_filters[j][k] = new LowPass (gamma_v[j],*lp_filters[j][k-1],x);
    }
  }
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "LPFilterGroup::~LPFilterGroup()"
LPFilterGroup::~LPFilterGroup() {
  for (int j=0; j<nalpha; j++) {
    for (int k=0; k<n_v[j]; k++) {
      delete lp_filters[j][k];
    }
  }
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void LPFilterGroup::print_some(const char *,Dofmap *,set<int> &)"
void LPFilterGroup::print_some(const char *filename,Dofmap *dofmap,
			       set<int> & node_list) {
  for (int j=0; j<nalpha; j++) {
    for (int k=0; k<n_v[j]; k++) {
      lp_filters[j][k]->state().print_some(filename,dofmap,node_list);
    }
  }
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void LPFilterGroup::update(Time time)"
void LPFilterGroup::update(Time time) {
  for (int j=0; j<nalpha; j++) {
    lp_filters[j][n_v[j]-1]->update(time);
    for (int k=0; k<n_v[j]; k++) {
      State &s = lp_filters[j][k]->state();
      s.set_time(time);
    }
  }
}  
