// -*-mode: c++ -*-
#ifndef SECANT_H
#define SECANT_H

#include <cassert>
#include <cmath>

/// Solves a 1D non-linear eq. with the secant method $f(x)=0$. Very simple. 
class Secant {
public:
  /// Initial value
  double x0;
  /// Initial increment 
  double epsilon;
  /// Tolerance for the solution of the equation
  double tol;
  /// function value at the solution
  double f;
  /// maximum number of iterations
  int maxits;
  /// iterations actually perfomed
  int its;
  /// relaxation factor
  double omega;
  /// Default constructor
  Secant() {};
  /// Constructor with default parameters
  Secant(double x0_) : x0(x0_), epsilon(1e-3), tol(1e-3), maxits(20),
    omega(1.) {};
  /// This should be defined by the user
  virtual double residual(double x) =0;
  /// Computes the solution
  double sol();
};
#endif
