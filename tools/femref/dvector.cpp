#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <libguile.h>

#include <src/dvector.h>

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
#define __FUN__ DVTYPE "-fun"
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
#define __FUN__ DVTYPE "-resize!"
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
#define __FUN__ DVTYPE "?"
static SCM
DVECTOR_P_FUN(SCM s_w) {
  return SCM_BOOL(SCM_SMOB_PREDICATE(TAG, s_w));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ DVTYPE "-set!"
static SCM
DVECTOR_SET_FUN(SCM s_w,SCM s_j,SCM s_v) {
  SCM_ASSERT (SCM_SMOB_PREDICATE(TAG, s_w),
              s_w, SCM_ARG1, __FUN__);
  dvector_t *w = (dvector_t *) SCM_SMOB_DATA(s_w);

  SCM_ASSERT(SCM_INUMP (s_j),s_j,SCM_ARG2,__FUN__);
  int j = SCM_INUM(s_j);
  
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

  w->ref(j) = v;
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

  SCM_ASSERT(SCM_INUMP (s_j),s_j,SCM_ARG2,__FUN__);
  int j = SCM_INUM(s_j);

  TYPE v = w->ref(j);
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
DVECTOR_DUMP_FUN(SCM s_w,SCM s_file) {
  dvector_t *w;

  SCM_ASSERT(SCM_SMOB_PREDICATE(TAG, s_w),
	     s_w, SCM_ARG1, __FUN__);
  w = (dvector_t *) SCM_SMOB_DATA (s_w);
  
  FILE *stream=NULL;
  if (s_file == SCM_UNDEFINED) {
    w->print();
  } else {
    SCM_ASSERT(scm_string_p(s_file),
	       s_file, SCM_ARG2, __FUN__);
    const char *file = SCM_STRING_CHARS(s_file);
    w->print(file);
  }
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

#define TT "#<" DVTYPE " "
  scm_puts (TT, port);
  int n = w->size();
  sprintf(buff,"%p, %d (",w,n);
  scm_lfwrite (buff,strlen(buff),port);
#define NPMAX 20
  int nn = (n>NPMAX ? NPMAX : n);
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
  if(n>NPMAX) 
    scm_puts (" ... ", port);
  scm_puts (")>", port);
  return 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
extern "C" void
INIT_DVECTOR_FUN(void) {
  TAG = scm_make_smob_type(DVTYPE, sizeof(dvector_t));
  scm_set_smob_free(TAG,FREE_DVECTOR_FUN);
  scm_set_smob_print (TAG,DVECTOR_PRINT_FUN);
  scm_c_define_gsubr("make-" DVTYPE, 0, 0, 0, scm_fun(MAKE_DVECTOR_FUN));
  scm_c_define_gsubr(DVTYPE "-push!", 2, 0, 0, scm_fun(DVECTOR_PUSH_FUN));
  scm_c_define_gsubr(DVTYPE "-size", 1, 0, 0, scm_fun(DVECTOR_SIZE_FUN));
  scm_c_define_gsubr(DVTYPE "-resize!", 2, 0, 0, scm_fun(DVECTOR_RESIZE_FUN));
  scm_c_define_gsubr(DVTYPE "-set!", 3, 0, 0, scm_fun(DVECTOR_SET_FUN));
  scm_c_define_gsubr(DVTYPE "-ref", 2, 0, 0, scm_fun(DVECTOR_REF_FUN));
  scm_c_define_gsubr(DVTYPE "-read!", 2, 0, 0, scm_fun(DVECTOR_READ_FUN));
  scm_c_define_gsubr(DVTYPE "-cat!", 2, 0, 0, scm_fun(DVECTOR_CAT_FUN));
  scm_c_define_gsubr(DVTYPE "-dump", 1, 1, 0, scm_fun(DVECTOR_DUMP_FUN));
}
