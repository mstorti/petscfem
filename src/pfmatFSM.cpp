#include "pfmatFSM.h"
static char _versID[] = "No Version.";
pfmatFSMin_assemblyState pfmatFSM::in_assemblyState;
pfmatFSMprofiledState pfmatFSM::profiledState;
pfmatFSMprofilingState pfmatFSM::profilingState;
pfmatFSMcleanState pfmatFSM::cleanState;
void pfmatFSMState::set_value(pfmatFSM& s)
  {s.FSMError("set_value", s.GetState().StateName());}
void pfmatFSMState::create(pfmatFSM& s)
  {s.FSMError("create", s.GetState().StateName());}
void pfmatFSMState::clear(pfmatFSM& s)
  {s.FSMError("clear", s.GetState().StateName());}
void pfmatFSMState::set_profile(pfmatFSM& s)
  {s.FSMError("set_profile", s.GetState().StateName());}
void pfmatFSMin_assemblyState::set_value(pfmatFSM& s) {
  s.SetState(pfmatFSM::in_assemblyState);
}
void pfmatFSMprofiledState::set_value(pfmatFSM& s) {
  s.SetState(pfmatFSM::in_assemblyState);
}
void pfmatFSMprofilingState::set_value(pfmatFSM& s) {
  s.SetState(pfmatFSM::profiledState);
  s.create();
  s.set_value();
}
void pfmatFSMprofilingState::create(pfmatFSM& s) {
  s.SetState(pfmatFSM::profiledState);
}
void pfmatFSMprofilingState::set_profile(pfmatFSM& s) {
  s.SetState(pfmatFSM::profilingState);
}
void pfmatFSMcleanState::clear(pfmatFSM& s) {
  s.SetState(pfmatFSM::cleanState);
}
void pfmatFSMcleanState::set_profile(pfmatFSM& s) {
  s.SetState(pfmatFSM::profilingState);
}
pfmatFSM::pfmatFSM() : itsState(&cleanState) {}
