#include "pfmatFSM.h"
static char _versID[] = "No Version.";
pfmatFSMfactoredState pfmatFSM::factoredState;
pfmatFSMassembledState pfmatFSM::assembledState;
pfmatFSMin_scatterState pfmatFSM::in_scatterState;
pfmatFSMin_assemblyState pfmatFSM::in_assemblyState;
pfmatFSMprofiledState pfmatFSM::profiledState;
pfmatFSMprofilingState pfmatFSM::profilingState;
pfmatFSMcleanState pfmatFSM::cleanState;
void pfmatFSMState::clean_factor(pfmatFSM& s)
  {s.FSMError("clean_factor", s.GetState().StateName());}
void pfmatFSMState::solve_only(pfmatFSM& s)
  {s.FSMError("solve_only", s.GetState().StateName());}
void pfmatFSMState::solve(pfmatFSM& s)
  {s.FSMError("solve", s.GetState().StateName());}
void pfmatFSMState::factor_and_solve(pfmatFSM& s)
  {s.FSMError("factor_and_solve", s.GetState().StateName());}
void pfmatFSMState::clean_mat(pfmatFSM& s)
  {s.FSMError("clean_mat", s.GetState().StateName());}
void pfmatFSMState::assembly_end(pfmatFSM& s)
  {s.FSMError("assembly_end", s.GetState().StateName());}
void pfmatFSMState::assembly_begin(pfmatFSM& s)
  {s.FSMError("assembly_begin", s.GetState().StateName());}
void pfmatFSMState::clean_prof(pfmatFSM& s)
  {s.FSMError("clean_prof", s.GetState().StateName());}
void pfmatFSMState::set_value(pfmatFSM& s)
  {s.FSMError("set_value", s.GetState().StateName());}
void pfmatFSMState::create(pfmatFSM& s)
  {s.FSMError("create", s.GetState().StateName());}
void pfmatFSMState::clear(pfmatFSM& s)
  {s.FSMError("clear", s.GetState().StateName());}
void pfmatFSMState::set_profile(pfmatFSM& s)
  {s.FSMError("set_profile", s.GetState().StateName());}
void pfmatFSMfactoredState::clear(pfmatFSM& s) {
  s.SetState(pfmatFSM::factoredState);
  s.clean_factor();
  s.clear();
}
void pfmatFSMfactoredState::clean_factor(pfmatFSM& s) {
  s.SetState(pfmatFSM::profiledState);
  s.clean_factor_a();
}
void pfmatFSMfactoredState::solve(pfmatFSM& s) {
  s.SetState(pfmatFSM::factoredState);
  s.solve_only_A();
}
void pfmatFSMfactoredState::solve_only(pfmatFSM& s) {
  s.SetState(pfmatFSM::factoredState);
  s.solve_only_A();
}
void pfmatFSMassembledState::solve(pfmatFSM& s) {
  s.SetState(pfmatFSM::factoredState);
  s.factor_and_solve_A();
}
void pfmatFSMassembledState::factor_and_solve(pfmatFSM& s) {
  s.SetState(pfmatFSM::factoredState);
  s.factor_and_solve_A();
}
void pfmatFSMassembledState::clear(pfmatFSM& s) {
  s.SetState(pfmatFSM::assembledState);
  s.clean_mat();
  s.clean_prof();
}
void pfmatFSMassembledState::clean_mat(pfmatFSM& s) {
  s.SetState(pfmatFSM::profiledState);
  s.clean_mat_a();
}
void pfmatFSMassembledState::set_value(pfmatFSM& s) {
  s.SetState(pfmatFSM::in_assemblyState);
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
void pfmatFSMprofiledState::set_value(pfmatFSM& s) {
  s.SetState(pfmatFSM::in_assemblyState);
}
void pfmatFSMprofiledState::clear(pfmatFSM& s) {
  s.SetState(pfmatFSM::profiledState);
  s.clean_mat();
}
void pfmatFSMprofiledState::clean_prof(pfmatFSM& s) {
  s.SetState(pfmatFSM::cleanState);
  s.clean_prof_a();
}
void pfmatFSMprofilingState::set_value(pfmatFSM& s) {
  s.SetState(pfmatFSM::profilingState);
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
