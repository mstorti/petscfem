// -*- mode: C++ -*-
//__INSERT_LICENSE__
//$Id: arglist.h,v 1.4 2001/04/01 01:35:06 mstorti Exp $

#ifndef ARGLIST_H
#define ARGLIST_H

#include <vector>
#include <string>
#include <mat.h>
#include <vec.h>
#include "libretto.h"

/// Type of arguments and processing information for items in the arglists's.
enum arg_options {

  // Operations to be performed
  DOWNLOAD_VECTOR        = 0x00000001,
  UPLOAD_VECTOR          = 0x00000002,
  ASSEMBLY_VECTOR        = 0x00000004,
  UPLOAD_MATRIX          = 0x00000008,
  ASSEMBLY_MATRIX        = 0x00000010,
  UPLOAD_PROFILE         = 0x00000020,
  ALLOC_MATRIX           = 0x00000040,
  UPLOAD_VECTOR_LOCST    = 0x00000080,
  IS_PERT_VECTOR         = 0x00000100,
  IS_FDJ_MATRIX          = 0x00000200,
  IS_FDJ_PROFILE         = 0x00000400,
  VECTOR_MIN             = 0x00000800,
  VECTOR_MAX             = 0x00001000,
  VECTOR_ADD             = 0x00002000,
  USER_DATA              = 0x00004000,

  // Some useful combinations. 
  IS_FDJ = IS_FDJ_PROFILE | IS_FDJ_MATRIX,
  UPLOAD_RETVAL = UPLOAD_VECTOR | UPLOAD_MATRIX
  | UPLOAD_PROFILE | UPLOAD_VECTOR_LOCST, 
  DELETE_RETVAL = UPLOAD_VECTOR | ALLOC_MATRIX| UPLOAD_PROFILE,
  VECTOR_ASSOC = VECTOR_MIN | VECTOR_MAX | VECTOR_ADD,

  // Type of args.
  IN_VECTOR = DOWNLOAD_VECTOR,
  PERT_VECTOR = DOWNLOAD_VECTOR | IS_PERT_VECTOR,
  IN_OUT_VECTOR = DOWNLOAD_VECTOR | UPLOAD_VECTOR_LOCST,
  OUT_VECTOR = UPLOAD_VECTOR | ASSEMBLY_VECTOR,
  OUT_MATRIX = ALLOC_MATRIX | UPLOAD_MATRIX | ASSEMBLY_MATRIX,
  OUT_MATRIX_FDJ = IS_FDJ_MATRIX | ASSEMBLY_MATRIX,
  PROFILE = ALLOC_MATRIX | UPLOAD_PROFILE,
  FDJ_PROFILE = IS_FDJ_PROFILE | UPLOAD_PROFILE
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Individual entry in the argument
    @author M. Storti
    @param arg (input) 
    @param arginfo (input) 
    
*/ 
class arg_entry {
public:
  /// Generic pointer to the argument being passed. 
  void *arg;

  /// other options in the form of bitfields. 
  int options;

  /** string containing info about what the
      argument contains and how has to be processed.
  */
  string arginfo;

  /// Constructor
  arg_entry(void *arg_,int options_, string arginfo_="") :
    arg(arg_), options(options_) {};
    
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/* @memo Arguments are passed to the `assemble' functions via 
    these `argument lists'. 
    @author M. Storti
    @param (input)
*/ 
class arg_list : public vector<arg_entry> {
public:
  void arg_add(void *arg,int options,string arginfo="");
};  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** For each argument the quantities passed to the element 
    routine are stored in a list of this structs.
    For instance for a MPI PETSc vector, we store a pointer 
    to it, a pointer to the `got' array of doubles in `sstate'
    and also a pointer to an array with ghost values. For an
    associative vector we store a pointer to the vector and a
    flag called 'was\_set' that flags if an element has been
    inspected or not. 
    @author M. Storti
*/ 
class arg_data {
public:
  /// A copy of the options for the corresponding arg\_entry value. 
  int options;
  /// The MPI vector.
  Vec *x;
  /// Vector of doubles where the MPI vector is `got'. 
  double *sstate;
  /// Sequential vector with ghost values. 
  Vec *ghost_vec;
  /** Vector of doubles where the ghost vector is `got' (with
      VecGetArray). 
  */
  double *ghost_vals;
  /// Local values (one row per element). 
  double *locst;
  /// MPI Matrix 
  Mat *A;
  /// Returned local values(one row per element) . 
  double *retval;
  /// reference state for finite difference calculation of jacobians.
  double *refres;
  /// Libretto dynamic array for storing the profile. 
  Darray *da;
  /// Vector for associative operationes (max, min, and add operations)
  vector<double> *vector_assoc;
  /** For vector\_assoc arguments with MIN and MAX  operations, flags
      if the actual min or max value has been already set. */
  int was_set;
  /** This is passed to the element routines and nothing specific is
      done. It is assumed that is managed by the user. */
  void *user_data;
  /// Default constructor. 
  arg_data(void);
};

typedef vector<arg_data> arg_data_list;

#endif

