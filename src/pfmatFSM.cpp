// Insert Debug Info for SMC Finite State Machine generated code
// ==insdebinfo== -- ./insdeb.pl
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
  printf("from: \"factored\", event: \"clear\","
         "to: \"factored\"\n");
  s.SetState(pfmatFSM::factoredState);
  s.clean_factor();
  s.clear();
}
void pfmatFSMfactoredState::clean_factor(pfmatFSM& s) {
  printf("from: \"factored\", event: \"clean_factor\","
         "to: \"assembled\"\n");
  s.SetState(pfmatFSM::assembledState);
  s.clean_factor_a();
  s.clean_mat_a();
}
void pfmatFSMfactoredState::solve(pfmatFSM& s) {
  printf("from: \"factored\", event: \"solve\","
         "to: \"factored\"\n");
  s.SetState(pfmatFSM::factoredState);
  s.solve_only_A();
}
void pfmatFSMfactoredState::solve_only(pfmatFSM& s) {
  printf("from: \"factored\", event: \"solve_only\","
         "to: \"factored\"\n");
  s.SetState(pfmatFSM::factoredState);
  s.solve_only_A();
}
void pfmatFSMassembledState::solve(pfmatFSM& s) {
  printf("from: \"assembled\", event: \"solve\","
         "to: \"factored\"\n");
  s.SetState(pfmatFSM::factoredState);
  s.factor_and_solve_A();
}
void pfmatFSMassembledState::factor_and_solve(pfmatFSM& s) {
  printf("from: \"assembled\", event: \"factor_and_solve\","
         "to: \"factored\"\n");
  s.SetState(pfmatFSM::factoredState);
  s.factor_and_solve_A();
}
void pfmatFSMassembledState::clear(pfmatFSM& s) {
  printf("from: \"assembled\", event: \"clear\","
         "to: \"assembled\"\n");
  s.SetState(pfmatFSM::assembledState);
  s.clean_mat();
  s.clean_prof();
}
void pfmatFSMassembledState::clean_mat(pfmatFSM& s) {
  printf("from: \"assembled\", event: \"clean_mat\","
         "to: \"profiled\"\n");
  s.SetState(pfmatFSM::profiledState);
  s.clean_mat_a();
}
void pfmatFSMassembledState::set_value(pfmatFSM& s) {
  printf("from: \"assembled\", event: \"set_value\","
         "to: \"in_assembly\"\n");
  s.SetState(pfmatFSM::in_assemblyState);
}
void pfmatFSMin_scatterState::assembly_end(pfmatFSM& s) {
  printf("from: \"in_scatter\", event: \"assembly_end\","
         "to: \"assembled\"\n");
  s.SetState(pfmatFSM::assembledState);
}
void pfmatFSMin_assemblyState::assembly_begin(pfmatFSM& s) {
  printf("from: \"in_assembly\", event: \"assembly_begin\","
         "to: \"in_scatter\"\n");
  s.SetState(pfmatFSM::in_scatterState);
}
void pfmatFSMin_assemblyState::set_value(pfmatFSM& s) {
  printf("from: \"in_assembly\", event: \"set_value\","
         "to: \"in_assembly\"\n");
  s.SetState(pfmatFSM::in_assemblyState);
}
void pfmatFSMprofiledState::set_value(pfmatFSM& s) {
  printf("from: \"profiled\", event: \"set_value\","
         "to: \"in_assembly\"\n");
  s.SetState(pfmatFSM::in_assemblyState);
}
void pfmatFSMprofiledState::clear(pfmatFSM& s) {
  printf("from: \"profiled\", event: \"clear\","
         "to: \"profiled\"\n");
  s.SetState(pfmatFSM::profiledState);
  s.clean_mat();
}
void pfmatFSMprofiledState::clean_prof(pfmatFSM& s) {
  printf("from: \"profiled\", event: \"clean_prof\","
         "to: \"clean\"\n");
  s.SetState(pfmatFSM::cleanState);
  s.clean_prof_a();
}
void pfmatFSMprofilingState::set_value(pfmatFSM& s) {
  printf("from: \"profiling\", event: \"set_value\","
         "to: \"profiling\"\n");
  s.SetState(pfmatFSM::profilingState);
  s.create();
  s.set_value();
}
void pfmatFSMprofilingState::create(pfmatFSM& s) {
  printf("from: \"profiling\", event: \"create\","
         "to: \"profiled\"\n");
  s.SetState(pfmatFSM::profiledState);
}
void pfmatFSMprofilingState::set_profile(pfmatFSM& s) {
  printf("from: \"profiling\", event: \"set_profile\","
         "to: \"profiling\"\n");
  s.SetState(pfmatFSM::profilingState);
}
void pfmatFSMcleanState::clear(pfmatFSM& s) {
  printf("from: \"clean\", event: \"clear\","
         "to: \"clean\"\n");
  s.SetState(pfmatFSM::cleanState);
}
void pfmatFSMcleanState::set_profile(pfmatFSM& s) {
  printf("from: \"clean\", event: \"set_profile\","
         "to: \"profiling\"\n");
  s.SetState(pfmatFSM::profilingState);
}
pfmatFSM::pfmatFSM() : itsState(&cleanState) {}
