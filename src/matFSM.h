#ifndef _H_MatFSM
#define _H_MatFSM
#include <stddef.h>
#include "sparse.h"
class MatFSM;

class MatFSMState {
public:

  virtual const char* StateName() const = 0;
  virtual void clear(MatFSM& s);
  virtual void solve(MatFSM& s);
  virtual void fill(MatFSM& s);
};

class MatFSMfactoredState : public MatFSMState {
public:
  virtual const char* StateName() const
  {return("factored");};
  virtual void clear(MatFSM&);
  virtual void solve(MatFSM&);
  virtual void fill(MatFSM&);
};

class MatFSMfilledState : public MatFSMState {
public:
  virtual const char* StateName() const
  {return("filled");};
  virtual void clear(MatFSM&);
  virtual void solve(MatFSM&);
  virtual void fill(MatFSM&);
};

class MatFSMcleanState : public MatFSMState {
public:
  virtual const char* StateName() const
  {return("clean");};
  virtual void clear(MatFSM&);
  virtual void fill(MatFSM&);
};
class MatFSM : public MatFSMContext {
  public:
  static MatFSMfactoredState factoredState;
  static MatFSMfilledState filledState;
  static MatFSMcleanState cleanState;
  MatFSM();// default constructor
  void clear() {itsState->clear(*this);}
  void solve() {itsState->solve(*this);}
  void fill() {itsState->fill(*this);}
  void SetState(MatFSMState& theState) {itsState=&theState;}
  MatFSMState& GetState() const {return *itsState;};
  private:
    MatFSMState* itsState;
};
#endif
