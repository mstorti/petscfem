#ifndef _H_pfmatFSM
#define _H_pfmatFSM
#include <stddef.h>
#include "pfmat.h"
class pfmatFSM;

class pfmatFSMState {
public:

  virtual const char* StateName() const = 0;
  virtual void set_value(pfmatFSM& s);
  virtual void create(pfmatFSM& s);
  virtual void clear(pfmatFSM& s);
  virtual void set_profile(pfmatFSM& s);
};

class pfmatFSMin_assemblyState : public pfmatFSMState {
public:
  virtual const char* StateName() const
  {return("in_assembly");};
  virtual void set_value(pfmatFSM&);
};

class pfmatFSMprofiledState : public pfmatFSMState {
public:
  virtual const char* StateName() const
  {return("profiled");};
  virtual void set_value(pfmatFSM&);
};

class pfmatFSMprofilingState : public pfmatFSMState {
public:
  virtual const char* StateName() const
  {return("profiling");};
  virtual void set_value(pfmatFSM&);
  virtual void create(pfmatFSM&);
  virtual void set_profile(pfmatFSM&);
};

class pfmatFSMcleanState : public pfmatFSMState {
public:
  virtual const char* StateName() const
  {return("clean");};
  virtual void clear(pfmatFSM&);
  virtual void set_profile(pfmatFSM&);
};
class pfmatFSM : public pfmatFSMContext {
  public:
  static pfmatFSMin_assemblyState in_assemblyState;
  static pfmatFSMprofiledState profiledState;
  static pfmatFSMprofilingState profilingState;
  static pfmatFSMcleanState cleanState;
  pfmatFSM();// default constructor
  void set_value() {itsState->set_value(*this);}
  void create() {itsState->create(*this);}
  void clear() {itsState->clear(*this);}
  void set_profile() {itsState->set_profile(*this);}
  void SetState(pfmatFSMState& theState) {itsState=&theState;}
  pfmatFSMState& GetState() const {return *itsState;};
  private:
    pfmatFSMState* itsState;
};
#endif
