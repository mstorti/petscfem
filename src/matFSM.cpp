#include "matFSM.h"
static char _versID[] = "No Version.";
MatFSMfactoredState MatFSM::factoredState;
MatFSMfilledState MatFSM::filledState;
MatFSMcleanState MatFSM::cleanState;
void MatFSMState::clear(MatFSM& s)
  {s.FSMError("clear", s.GetState().StateName());}
void MatFSMState::solve(MatFSM& s)
  {s.FSMError("solve", s.GetState().StateName());}
void MatFSMState::fill(MatFSM& s)
  {s.FSMError("fill", s.GetState().StateName());}
void MatFSMfactoredState::clear(MatFSM& s) {
  s.SetState(MatFSM::filledState);
  s.clean_factor();
  s.clear();
}
void MatFSMfactoredState::solve(MatFSM& s) {
  s.SetState(MatFSM::factoredState);
  s.solve_only();
}
void MatFSMfactoredState::fill(MatFSM& s) {
  s.SetState(MatFSM::filledState);
  s.clean_factor();
}
void MatFSMfilledState::clear(MatFSM& s) {
  s.SetState(MatFSM::cleanState);
  s.clean_mat();
}
void MatFSMfilledState::solve(MatFSM& s) {
  s.SetState(MatFSM::factoredState);
  s.fact_and_solve();
}
void MatFSMfilledState::fill(MatFSM& s) {
  s.SetState(MatFSM::filledState);
}
void MatFSMcleanState::clear(MatFSM& s) {
  s.SetState(MatFSM::cleanState);
}
void MatFSMcleanState::fill(MatFSM& s) {
  s.SetState(MatFSM::filledState);
}
MatFSM::MatFSM() : itsState(&cleanState) {}
