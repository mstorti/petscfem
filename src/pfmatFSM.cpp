#include "pfmatFSM.h"
static char _versID[] = "No Version.";
pfmatFSMfactoredState pfmatFSM::factoredState;
pfmatFSMassembledState pfmatFSM::assembledState;
pfmatFSMin_scatterState pfmatFSM::in_scatterState;
pfmatFSMin_assemblyState pfmatFSM::in_assemblyState;
pfmatFSMprofiledState pfmatFSM::profiledState;
pfmatFSMcleanState pfmatFSM::cleanState;
void pfmatFSMState::solve(pfmatFSM& s)
  {s.FSMError("solve", s.GetState().StateName());}
void pfmatFSMState::zero_entries(pfmatFSM& s)
  {s.FSMError("zero_entries", s.GetState().StateName());}
void pfmatFSMState::assembly_end(pfmatFSM& s)
  {s.FSMError("assembly_end", s.GetState().StateName());}
void pfmatFSMState::assembly_begin(pfmatFSM& s)
  {s.FSMError("assembly_begin", s.GetState().StateName());}
void pfmatFSMState::clear(pfmatFSM& s)
  {s.FSMError("clear", s.GetState().StateName());}
void pfmatFSMState::set_value(pfmatFSM& s)
  {s.FSMError("set_value", s.GetState().StateName());}
void pfmatFSMState::set_profile(pfmatFSM& s)
  {s.FSMError("set_profile", s.GetState().StateName());}
void pfmatFSMfactoredState::clear(pfmatFSM& s) {
  s.SetState(pfmatFSM::cleanState);
  s.clean_factor();
}
void pfmatFSMfactoredState::zero_entries(pfmatFSM& s) {
  s.SetState(pfmatFSM::profiledState);
  s.clean_factor();
}
void pfmatFSMfactoredState::solve(pfmatFSM& s) {
  s.SetState(pfmatFSM::factoredState);
  s.solve_only();
}
void pfmatFSMassembledState::solve(pfmatFSM& s) {
  s.SetState(pfmatFSM::factoredState);
  s.factor_and_solve();
}
void pfmatFSMassembledState::clear(pfmatFSM& s) {
  s.SetState(pfmatFSM::cleanState);
}
void pfmatFSMassembledState::set_value(pfmatFSM& s) {
  s.SetState(pfmatFSM::in_assemblyState);
}
void pfmatFSMassembledState::zero_entries(pfmatFSM& s) {
  s.SetState(pfmatFSM::profiledState);
}
void pfmatFSMin_scatterState::assembly_end(pfmatFSM& s) {
  s.SetState(pfmatFSM::assembledState);
}
void pfmatFSMin_assemblyState::assembly_begin(pfmatFSM& s) {
  s.SetState(pfmatFSM::in_scatterState);
}
void pfmatFSMin_assemblyState::set_value(pfmatFSM& s) {
  s.SetState(pfmatFSM::in_assemblyState);
}
void pfmatFSMprofiledState::clear(pfmatFSM& s) {
  s.SetState(pfmatFSM::cleanState);
  s.clean_profile();
}
void pfmatFSMprofiledState::set_value(pfmatFSM& s) {
  s.SetState(pfmatFSM::in_assemblyState);
}
void pfmatFSMcleanState::set_profile(pfmatFSM& s) {
  s.SetState(pfmatFSM::profiledState);
}
pfmatFSM::pfmatFSM() : itsState(&cleanState) {}
