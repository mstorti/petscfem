#ifndef _H_pfmatFSM
#define _H_pfmatFSM
#include <stddef.h>
#include "pfmat.h"
class pfmatFSM;

class pfmatFSMState {
public:

  virtual const char* StateName() const = 0;
  virtual void solve(pfmatFSM& s);
  virtual void zero_entries(pfmatFSM& s);
  virtual void assembly_end(pfmatFSM& s);
  virtual void assembly_begin(pfmatFSM& s);
  virtual void clear(pfmatFSM& s);
  virtual void set_value(pfmatFSM& s);
  virtual void set_profile(pfmatFSM& s);
};

class pfmatFSMfactoredState : public pfmatFSMState {
public:
  virtual const char* StateName() const
  {return("factored");};
  virtual void clear(pfmatFSM&);
  virtual void zero_entries(pfmatFSM&);
  virtual void solve(pfmatFSM&);
};

class pfmatFSMassembledState : public pfmatFSMState {
public:
  virtual const char* StateName() const
  {return("assembled");};
  virtual void solve(pfmatFSM&);
  virtual void clear(pfmatFSM&);
  virtual void set_value(pfmatFSM&);
  virtual void zero_entries(pfmatFSM&);
};

class pfmatFSMin_scatterState : public pfmatFSMState {
public:
  virtual const char* StateName() const
  {return("in_scatter");};
  virtual void assembly_end(pfmatFSM&);
};

class pfmatFSMin_assemblyState : public pfmatFSMState {
public:
  virtual const char* StateName() const
  {return("in_assembly");};
  virtual void assembly_begin(pfmatFSM&);
  virtual void set_value(pfmatFSM&);
};

class pfmatFSMprofiledState : public pfmatFSMState {
public:
  virtual const char* StateName() const
  {return("profiled");};
  virtual void clear(pfmatFSM&);
  virtual void set_value(pfmatFSM&);
};

class pfmatFSMcleanState : public pfmatFSMState {
public:
  virtual const char* StateName() const
  {return("clean");};
  virtual void set_profile(pfmatFSM&);
};
class pfmatFSM : public pfmatFSMContext {
  public:
  static pfmatFSMfactoredState factoredState;
  static pfmatFSMassembledState assembledState;
  static pfmatFSMin_scatterState in_scatterState;
  static pfmatFSMin_assemblyState in_assemblyState;
  static pfmatFSMprofiledState profiledState;
  static pfmatFSMcleanState cleanState;
  pfmatFSM();// default constructor
  void solve() {itsState->solve(*this);}
  void zero_entries() {itsState->zero_entries(*this);}
  void assembly_end() {itsState->assembly_end(*this);}
  void assembly_begin() {itsState->assembly_begin(*this);}
  void clear() {itsState->clear(*this);}
  void set_value() {itsState->set_value(*this);}
  void set_profile() {itsState->set_profile(*this);}
  void SetState(pfmatFSMState& theState) {itsState=&theState;}
  pfmatFSMState& GetState() const {return *itsState;};
  private:
    pfmatFSMState* itsState;
};
#endif
