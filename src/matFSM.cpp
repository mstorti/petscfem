#include "matFSM.h"
static char _versID[] = "No Version.";
MatFSMfactoredState MatFSM::factoredState;
MatFSMfilledState MatFSM::filledState;
MatFSMcleanState MatFSM::cleanState;
void MatFSMState::clear(MatFSM& s)
  {s.FSMError("clear", s.GetState().StateName());}
void MatFSMState::factor(MatFSM& s)
  {s.FSMError("factor", s.GetState().StateName());}
void MatFSMState::solve(MatFSM& s)
  {s.FSMError("solve", s.GetState().StateName());}
void MatFSMState::fill(MatFSM& s)
  {s.FSMError("fill", s.GetState().StateName());}
void MatFSMfactoredState::clear(MatFSM& s) {
  s.SetState(MatFSM::cleanState);
  s.clean_factor();
  s.clean_mat();
}
void MatFSMfactoredState::fill(MatFSM& s) {
  s.SetState(MatFSM::filledState);
  s.clean_factor();
}
void MatFSMfactoredState::solve(MatFSM& s) {
  s.SetState(MatFSM::factoredState);
  s.back_subst();
}
void MatFSMfilledState::clear(MatFSM& s) {
  s.SetState(MatFSM::cleanState);
  s.clean_mat();
}
void MatFSMfilledState::factor(MatFSM& s) {
  s.SetState(MatFSM::factoredState);
  s.factor();
}
void MatFSMfilledState::solve(MatFSM& s) {
  s.SetState(MatFSM::factoredState);
  s.factor();
  s.back_subst();
}
void MatFSMcleanState::fill(MatFSM& s) {
  s.SetState(MatFSM::filledState);
}
MatFSM::MatFSM() : itsState(&cleanState) {}
