// Insert Debug Info for SMC Finite State Machine generated code
// script version:  %Id: insdeb.pl,v 1.2 2002/01/14 03:45:06 mstorti Exp $ 
#include "pfmatFSM.h"
static char _versID[] = "No Version.";
pfmatFSMfactoredState pfmatFSM::factoredState;
pfmatFSMassembledState pfmatFSM::assembledState;
pfmatFSMin_assemblyState pfmatFSM::in_assemblyState;
pfmatFSMin_scatterState pfmatFSM::in_scatterState;
pfmatFSMprofiledState pfmatFSM::profiledState;
pfmatFSMprofilingState pfmatFSM::profilingState;
pfmatFSMcleanState pfmatFSM::cleanState;
void pfmatFSMState::solve_only(pfmatFSM& s)
  {s.FSMError("solve_only", s.GetState().StateName());}
void pfmatFSMState::solve(pfmatFSM& s)
  {s.FSMError("solve", s.GetState().StateName());}
void pfmatFSMState::factor_and_solve(pfmatFSM& s)
  {s.FSMError("factor_and_solve", s.GetState().StateName());}
void pfmatFSMState::assembly_end(pfmatFSM& s)
  {s.FSMError("assembly_end", s.GetState().StateName());}
void pfmatFSMState::assembly_begin(pfmatFSM& s)
  {s.FSMError("assembly_begin", s.GetState().StateName());}
void pfmatFSMState::clean_prof(pfmatFSM& s)
  {s.FSMError("clean_prof", s.GetState().StateName());}
void pfmatFSMState::clean_mat(pfmatFSM& s)
  {s.FSMError("clean_mat", s.GetState().StateName());}
void pfmatFSMState::clean_factor(pfmatFSM& s)
  {s.FSMError("clean_factor", s.GetState().StateName());}
void pfmatFSMState::set_value(pfmatFSM& s)
  {s.FSMError("set_value", s.GetState().StateName());}
void pfmatFSMState::clear(pfmatFSM& s)
  {s.FSMError("clear", s.GetState().StateName());}
void pfmatFSMState::asssembly_begin(pfmatFSM& s)
  {s.FSMError("asssembly_begin", s.GetState().StateName());}
void pfmatFSMState::create(pfmatFSM& s)
  {s.FSMError("create", s.GetState().StateName());}
void pfmatFSMState::set_profile(pfmatFSM& s)
  {s.FSMError("set_profile", s.GetState().StateName());}
void pfmatFSMfactoredState::clear(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("factored","clear","factored");
  s.SetState(pfmatFSM::factoredState);
  s.clean_factor();
  s.clear();
}
void pfmatFSMfactoredState::set_value(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("factored","set_value","factored");
  s.SetState(pfmatFSM::factoredState);
  s.clean_factor();
}
void pfmatFSMfactoredState::clean_factor(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("factored","clean_factor","assembled");
  s.SetState(pfmatFSM::assembledState);
  s.clean_factor_a();
  s.clean_mat_a();
}
void pfmatFSMfactoredState::solve(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("factored","solve","factored");
  s.SetState(pfmatFSM::factoredState);
  s.solve_only_A();
}
void pfmatFSMfactoredState::solve_only(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("factored","solve_only","factored");
  s.SetState(pfmatFSM::factoredState);
  s.solve_only_A();
}
void pfmatFSMfactoredState::clean_mat(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("factored","clean_mat","factored");
  s.SetState(pfmatFSM::factoredState);
  s.clean_factor();
}
void pfmatFSMassembledState::assembly_begin(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("assembled","assembly_begin","in_scatter");
  s.SetState(pfmatFSM::in_scatterState);
}
void pfmatFSMassembledState::solve(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("assembled","solve","factored");
  s.SetState(pfmatFSM::factoredState);
  s.factor_and_solve_A();
}
void pfmatFSMassembledState::factor_and_solve(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("assembled","factor_and_solve","factored");
  s.SetState(pfmatFSM::factoredState);
  s.factor_and_solve_A();
}
void pfmatFSMassembledState::clear(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("assembled","clear","assembled");
  s.SetState(pfmatFSM::assembledState);
  s.clean_mat();
  s.clean_prof();
}
void pfmatFSMassembledState::clean_mat(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("assembled","clean_mat","profiled");
  s.SetState(pfmatFSM::profiledState);
  s.clean_mat_a();
}
void pfmatFSMassembledState::set_value(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("assembled","set_value","in_assembly");
  s.SetState(pfmatFSM::in_assemblyState);
}
void pfmatFSMin_scatterState::assembly_end(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("in_scatter","assembly_end","assembled");
  s.SetState(pfmatFSM::assembledState);
}
void pfmatFSMin_assemblyState::clean_mat(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("in_assembly","clean_mat","in_assembly");
  s.SetState(pfmatFSM::in_assemblyState);
  s.assembly_begin();
  s.assembly_end();
  s.clean_mat();
}
void pfmatFSMin_assemblyState::assembly_begin(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("in_assembly","assembly_begin","in_scatter");
  s.SetState(pfmatFSM::in_scatterState);
}
void pfmatFSMin_assemblyState::set_value(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("in_assembly","set_value","in_assembly");
  s.SetState(pfmatFSM::in_assemblyState);
}
void pfmatFSMprofiledState::assembly_begin(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("profiled","assembly_begin","in_scatter");
  s.SetState(pfmatFSM::in_scatterState);
}
void pfmatFSMprofiledState::set_value(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("profiled","set_value","in_assembly");
  s.SetState(pfmatFSM::in_assemblyState);
}
void pfmatFSMprofiledState::clear(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("profiled","clear","profiled");
  s.SetState(pfmatFSM::profiledState);
  s.clean_mat();
}
void pfmatFSMprofiledState::clean_prof(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("profiled","clean_prof","clean");
  s.SetState(pfmatFSM::cleanState);
  s.clean_prof_a();
}
void pfmatFSMprofiledState::clean_mat(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("profiled","clean_mat","profiled");
  s.SetState(pfmatFSM::profiledState);
  s.clean_mat_a();
}
void pfmatFSMprofiledState::clean_factor(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("profiled","clean_factor","profiled");
  s.SetState(pfmatFSM::profiledState);
}
void pfmatFSMprofilingState::set_value(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("profiling","set_value","profiling");
  s.SetState(pfmatFSM::profilingState);
  s.create();
  s.set_value();
}
void pfmatFSMprofilingState::create(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("profiling","create","profiled");
  s.SetState(pfmatFSM::profiledState);
}
void pfmatFSMprofilingState::set_profile(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("profiling","set_profile","profiling");
  s.SetState(pfmatFSM::profilingState);
}
void pfmatFSMcleanState::clear(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("clean","clear","clean");
  s.SetState(pfmatFSM::cleanState);
}
void pfmatFSMcleanState::asssembly_begin(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("clean","asssembly_begin","in_scatter");
  s.SetState(pfmatFSM::in_scatterState);
}
void pfmatFSMcleanState::create(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("clean","create","profiled");
  s.SetState(pfmatFSM::profiledState);
}
void pfmatFSMcleanState::set_profile(pfmatFSM& s) {
  PRINT_FSM_TRANSITION_INFO("clean","set_profile","profiling");
  s.SetState(pfmatFSM::profilingState);
}
pfmatFSM::pfmatFSM() : itsState(&cleanState) {}
