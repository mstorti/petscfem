#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <libguile.h>

#include <src/dvector.h>
#include "./dvector2.h"
#include "./guilemac.h"

scm_t_bits TAG;

typedef SCM(*scm_fun)();

typedef dvector<TYPE> dvector_t;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
static SCM
MAKE_DVECTOR_FUN() {
  dvector_t *w;
  w = new dvector_t;
  SCM_RETURN_NEWSMOB (TAG, w);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ DVTYPE "(free)"
static size_t
FREE_DVECTOR_FUN(SCM s_w) {
  SCM_ASSERT (SCM_SMOB_PREDICATE(TAG, s_w),
              s_w, SCM_ARG1, __FUN__);
  dvector_t *w = (dvector_t *) SCM_SMOB_DATA(s_w);
  size_t size;
  size = 0;
  w->clear();
  delete w;
  return size;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ DVTYPE "-clone!"
static SCM
DVECTOR_CLONE_FUN(SCM s_v,SCM s_w) {

  SCM_ASSERT (SCM_SMOB_PREDICATE(TAG, s_v),
              s_v, SCM_ARG1, __FUN__);
  dvector_t *v = (dvector_t *) SCM_SMOB_DATA(s_v);

  SCM_ASSERT (SCM_SMOB_PREDICATE(TAG, s_w),
              s_w, SCM_ARG1, __FUN__);
  dvector_t *w = (dvector_t *) SCM_SMOB_DATA(s_w);

  v->clone(*w);
  return SCM_UNSPECIFIED;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if ! defined DV_INT
SCM_DEFINE(DVECTOR_SCALE_FUN,DVECTOR_SCALE_FUN_S, 2, 0, 0,
	   (SCM s_w,
	    SCM s_alpha),
	   "Scales w by alpha.")
#define FUNC_NAME S(DVECTOR_SCALE_FUN)
{
  DVDBLARG(w,1);
  MY_SCM_GET_DBL(alpha,2,(DVECTOR_SCALE_FUN_S ". Getting alpha. "));
  w->scale(alpha);
  return SCM_UNSPECIFIED;
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ DVTYPE "-resize-w!"
static SCM
DVECTOR_RESIZE_FUN(SCM s_w,SCM s_n) {
  SCM_ASSERT (SCM_SMOB_PREDICATE(TAG, s_w),
              s_w, SCM_ARG1, __FUN__);
  dvector_t *w = (dvector_t *) SCM_SMOB_DATA(s_w);

  SCM_ASSERT(SCM_INUMP (s_n),s_n,SCM_ARG2,__FUN__);
  int n = SCM_INUM(s_n);

  w->resize(n);
  return SCM_UNSPECIFIED;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ DVTYPE "-reshape!"
static SCM
DVECTOR_RESHAPE_FUN(SCM s_w,SCM s_shape) {
  SCM_ASSERT (SCM_SMOB_PREDICATE(TAG, s_w),
              s_w, SCM_ARG1, __FUN__);
  dvector_t *w = (dvector_t *) SCM_SMOB_DATA(s_w);
  vector<int> shape;
  scmlist2vec(s_shape,shape);
  w->reshape(shape);
  return SCM_UNSPECIFIED;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ DVTYPE "-shape"
static SCM
DVECTOR_SHAPE_FUN(SCM s_w) {
  SCM_ASSERT (SCM_SMOB_PREDICATE(TAG, s_w),
              s_w, SCM_ARG1, __FUN__);
  dvector_t *w = (dvector_t *) SCM_SMOB_DATA(s_w);
  vector<int> shape;
  int n = w->rank();
  shape.resize(n);
  for (int j=0; j<n; j++)
    shape[j] = w->size(j);
  return vec2scmlist(shape);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ DVTYPE "?"
static SCM
DVECTOR_P_FUN(SCM s_w) {
  return SCM_BOOL(SCM_SMOB_PREDICATE(TAG, s_w));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ DVTYPE "-set-w1"
static SCM
DVECTOR_SET_W1_FUN(SCM s_w,SCM s_v) {
  SCM_ASSERT (SCM_SMOB_PREDICATE(TAG, s_w),
              s_w, SCM_ARG1, __FUN__);
  dvector_t *w = (dvector_t *) SCM_SMOB_DATA(s_w);

  TYPE v;
#if defined DV_INT
  SCM_ASSERT(SCM_INUMP (s_v),s_v,SCM_ARG3,__FUN__);
  v = SCM_INUM(s_v);
#elif defined DV_DBL
  SCM_ASSERT(scm_real_p (s_v),s_v,SCM_ARG3, __FUN__);
  v = scm_num2dbl(s_v,__FUN__);
#else
#error undefined type!! 
#endif

  w->set(v);
  return SCM_UNSPECIFIED;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ DVTYPE "-set-w2"
static SCM
DVECTOR_SET_W2_FUN(SCM s_w,SCM s_j,SCM s_v) {
  SCM_ASSERT (SCM_SMOB_PREDICATE(TAG, s_w),
              s_w, SCM_ARG1, __FUN__);
  dvector_t *w = (dvector_t *) SCM_SMOB_DATA(s_w);

  TYPE v;
#if defined DV_INT
  SCM_ASSERT(SCM_INUMP (s_v),s_v,SCM_ARG3,__FUN__);
  v = SCM_INUM(s_v);
#elif defined DV_DBL
  SCM_ASSERT(scm_real_p (s_v),s_v,SCM_ARG3, __FUN__);
  v = scm_num2dbl(s_v,__FUN__);
#else
#error undefined type!! 
#endif

  if (SCM_INUMP(s_j)) {
    int j = SCM_INUM(s_j);
    w->ref(j) = v;
  } else if (scm_list_p(s_j)) {
    vector<int> indx;
    scmlist2vec(s_j,indx);
    w->e(indx) = v;
  } else {
    scm_wrong_type_arg(__FUN__,2,s_j);
  }
  return SCM_UNSPECIFIED;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ DVTYPE "-ref"
static SCM
DVECTOR_REF_FUN(SCM s_w,SCM s_j) {
  SCM_ASSERT (SCM_SMOB_PREDICATE(TAG, s_w),
              s_w, SCM_ARG1, __FUN__);
  dvector_t *w = (dvector_t *) SCM_SMOB_DATA(s_w);

  TYPE v;
  if (SCM_INUMP(s_j)) {
    int j = SCM_INUM(s_j);
    v = w->ref(j);
  } else if (scm_list_p(s_j)) {
    vector<int> indx;
    scmlist2vec(s_j,indx);
    v = w->e(indx);
  } else {
    scm_wrong_type_arg(__FUN__,2,s_j);
  }

#if defined DV_INT
  return SCM_MAKINUM(v);
#elif defined DV_DBL
  return scm_make_real(v);
#else
#error undefined type!! 
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ DVTYPE "-push!"
static SCM
DVECTOR_PUSH_FUN(SCM s_w, SCM s_v) {
  dvector_t *w;
  TYPE v;

  SCM_ASSERT (SCM_SMOB_PREDICATE(TAG, s_w),
              s_w, SCM_ARG1, __FUN__);
  w = (dvector_t *) SCM_SMOB_DATA (s_w);

#if defined DV_INT
  SCM_ASSERT(SCM_INUMP (s_v),s_v,SCM_ARG2,__FUN__);
  v = SCM_INUM(s_v);
#elif defined DV_DBL
  SCM_ASSERT(scm_real_p (s_v),s_v,SCM_ARG2, __FUN__);
  v = scm_num2dbl(s_v,__FUN__);
#else
#error undefined type!! 
#endif

  w->push(v);
  return SCM_UNSPECIFIED;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ DVTYPE "-size"
static SCM
DVECTOR_SIZE_FUN(SCM s_w) {
  dvector_t *w;

  SCM_ASSERT (SCM_SMOB_PREDICATE(TAG, s_w),
              s_w, SCM_ARG1, __FUN__);
  w = (dvector_t *) SCM_SMOB_DATA (s_w);

  return SCM_MAKINUM(w->size());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ DVTYPE "-read!"
static SCM
DVECTOR_READ_FUN(SCM s_w,SCM s_file) {
  dvector_t *w;

  SCM_ASSERT(SCM_SMOB_PREDICATE(TAG, s_w),
	     s_w, SCM_ARG1, __FUN__);
  w = (dvector_t *) SCM_SMOB_DATA (s_w);
  
  SCM_ASSERT(scm_string_p(s_file),
	     s_file, SCM_ARG2, __FUN__);
  
  const char *file = SCM_STRING_CHARS(s_file);
  w->read(file);
  return SCM_MAKINUM(w->size());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ DVTYPE "-cat!"
static SCM
DVECTOR_CAT_FUN(SCM s_w,SCM s_file) {
  dvector_t *w;

  SCM_ASSERT(SCM_SMOB_PREDICATE(TAG, s_w),
	     s_w, SCM_ARG1, __FUN__);
  w = (dvector_t *) SCM_SMOB_DATA (s_w);
  
  SCM_ASSERT(scm_string_p(s_file),
	     s_file, SCM_ARG2, __FUN__);
  
  const char *file = SCM_STRING_CHARS(s_file);
  w->cat(file);
  return SCM_MAKINUM(w->size());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ DVTYPE "-dump"
static SCM
DVECTOR_DUMP_FUN(SCM s_w,SCM s_file,SCM s_rowsz) {
  dvector_t *w;

  SCM_ASSERT(SCM_SMOB_PREDICATE(TAG, s_w),
	     s_w, SCM_ARG1, __FUN__);
  w = (dvector_t *) SCM_SMOB_DATA (s_w);

  int rowsz;
  if (s_rowsz == SCM_UNDEFINED) 
    rowsz = 0;
  else {
    SCM_ASSERT(SCM_INUMP(s_rowsz),
	       s_rowsz, SCM_ARG3, __FUN__);
    rowsz = SCM_INUM(s_rowsz);
  }

  bool tostdout = false;
  if (s_file == SCM_UNDEFINED) 
    tostdout = true;
  else {
    SCM_ASSERT(scm_string_p(s_file),
	       s_file, SCM_ARG2, __FUN__);
    if (scm_string_null_p(s_file))
      tostdout = true;
    else {
      const char *file = SCM_STRING_CHARS(s_file);
      w->print(file,rowsz);
    }
  }
  if (tostdout) w->print(stdout,rowsz);

  return SCM_UNSPECIFIED;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ DVTYPE "-push!"
static int
DVECTOR_PRINT_FUN(SCM s_w,SCM port, scm_print_state *pstate) {
  static char buff[100];
  dvector_t *w;
  int k;

  SCM_ASSERT (SCM_SMOB_PREDICATE (TAG, s_w),
              s_w, SCM_ARG1, __FUN__);
  w = (dvector_t *) SCM_SMOB_DATA (s_w);

#define TT "#,(" DVTYPE " "
  scm_puts (TT, port);
  int n = w->size();
  sprintf(buff," %d (",n);
  scm_lfwrite (buff,strlen(buff),port);
  int m = w->rank();
  for (int j=0; j<m; j++) {
    sprintf(buff,"%d ",w->size(j));
    scm_lfwrite (buff,strlen(buff),port);
  }
  sprintf(buff,") ");
  scm_lfwrite (buff,strlen(buff),port);
#define NPMAX 20
  // int nn = (n>NPMAX ? NPMAX : n);
  int nn = n;
  for (int j=0; j<nn; j++) {
#if defined DV_INT
#define PRINTF_FORMAT "%d "
#elif defined DV_DBL
#define PRINTF_FORMAT "%.12g "
#else
#error undefined type!! 
#endif
    sprintf(buff,PRINTF_FORMAT,w->ref(j));
    scm_lfwrite (buff,strlen(buff),port);
  }
  if(n>nn) 
    scm_puts (" ... ", port);
  scm_puts (")", port);
  return 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ DVTYPE "-read-ctor"
static int
DVECTOR_READ_CTOR_FUN(SCM port) {
  assert(0);
#if 0
  static char buff[100];
  dvector_t *w;
  int k;

  SCM_ASSERT (SCM_SMOB_PREDICATE (TAG, s_w),
              s_w, SCM_ARG1, __FUN__);
  w = (dvector_t *) SCM_SMOB_DATA (s_w);

#define TT "#,(" DVTYPE " "
  scm_puts (TT, port);
  int n = w->size();
  sprintf(buff," %d (",n);
  scm_lfwrite (buff,strlen(buff),port);
  int m = w->rank();
  for (int j=0; j<m; j++) {
    sprintf(buff,"%d ",w->size(j));
    scm_lfwrite (buff,strlen(buff),port);
  }
  sprintf(buff,") ");
  scm_lfwrite (buff,strlen(buff),port);
#define NPMAX 20
  // int nn = (n>NPMAX ? NPMAX : n);
  int nn = n;
  for (int j=0; j<nn; j++) {
#if defined DV_INT
#define PRINTF_FORMAT "%d "
#elif defined DV_DBL
#define PRINTF_FORMAT "%.12g "
#else
#error undefined type!! 
#endif
    sprintf(buff,PRINTF_FORMAT,w->ref(j));
    scm_lfwrite (buff,strlen(buff),port);
  }
  if(n>nn) 
    scm_puts (" ... ", port);
  scm_puts (")", port);
#endif
  return 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
extern "C" void
INIT_DVECTOR_FUN(void) {
  TAG = scm_make_smob_type(DVTYPE, sizeof(dvector_t));
  scm_set_smob_free(TAG,FREE_DVECTOR_FUN);
  scm_set_smob_print (TAG,DVECTOR_PRINT_FUN);
  scm_c_define_gsubr("make-" DVTYPE, 0, 0, 0, scm_fun(MAKE_DVECTOR_FUN));
  scm_c_define_gsubr(DVTYPE "-clone!", 2, 0, 0, scm_fun(DVECTOR_CLONE_FUN));
  scm_c_define_gsubr(DVTYPE "-push!", 2, 0, 0, scm_fun(DVECTOR_PUSH_FUN));
  scm_c_define_gsubr(DVTYPE "-size", 1, 0, 0, scm_fun(DVECTOR_SIZE_FUN));
  scm_c_define_gsubr(DVTYPE "-resize-w!", 2, 0, 0, scm_fun(DVECTOR_RESIZE_FUN));
  scm_c_define_gsubr(DVTYPE "-reshape!", 1, 0, 1, scm_fun(DVECTOR_RESHAPE_FUN));
  scm_c_define_gsubr(DVTYPE "-shape", 1, 0, 0, scm_fun(DVECTOR_SHAPE_FUN));
  scm_c_define_gsubr(DVTYPE "-set-w1", 2, 0, 0, scm_fun(DVECTOR_SET_W1_FUN));
  scm_c_define_gsubr(DVTYPE "-set-w2", 3, 0, 0, scm_fun(DVECTOR_SET_W2_FUN));
  scm_c_define_gsubr(DVTYPE "-ref", 2, 0, 0, scm_fun(DVECTOR_REF_FUN));
  scm_c_define_gsubr(DVTYPE "-read!", 2, 0, 0, scm_fun(DVECTOR_READ_FUN));
  scm_c_define_gsubr(DVTYPE "-cat!", 2, 0, 0, scm_fun(DVECTOR_CAT_FUN));
  scm_c_define_gsubr(DVTYPE "-dump", 1, 2, 0, scm_fun(DVECTOR_DUMP_FUN));
#ifndef SCM_MAGIC_SNARFER
#if defined DV_DBL
#include "./dvectord.x"
#elif defined DV_INT
#include "./dvectori.x"
#else
#error "not defined dv-type"
#endif
#endif

}
