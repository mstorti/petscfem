//__INSERT_LICENSE__
//$Id$

#include <src/debug.h>
#include <set>

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/pfmat.h>
#include <src/hook.h>

#include "advective.h"

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 

int  print_internal_loop_conv_g         = 0;
int  consistent_supg_matrix_g           = 0;
int  local_time_step_g                  = 0;
int  comp_mat_each_time_step_g          = 0;
int  verify_jacobian_with_numerical_one = 0;

const char * jobinfo_fields = 0;

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 

// Low-Mach stratgy
double local_dt_min;
