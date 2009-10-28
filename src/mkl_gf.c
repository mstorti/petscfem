/*
 * MKL-11 only? _gfortran_internal_malloc64 problem.
 *
 * MKL_LIB_PATH=/usr/local/mkl
 * ORIG_MKL_LIB_PATH=/opt/intel/Compiler/11.0/ * /mkl/lib/em64t/

 * sudo mkdir $MKL_LIB_PATH
 * sudo cp $ORIG_MKL_LIB_PATH/{libmkl_gf_lp64,libmkl_core}.a $MKL_LIB_PATH
 * wget http://prs.ism.ac.jp/~nakama/mkl/mkl_gf.c
 * gcc -O3 -fPIC -g -c mkl_gf.c
 * sudo cp $ORIG_MKL_LIB_PATH/libmkl_gnu_thread.a $MKL_LIB_PATH
 * sudo ar r $MKL_LIB_PATH/libmkl_gnu_thread.a mkl_gf.o
 *
 * MKL="   -L${MKL_LIB_PATH}                              \
 *         -Wl,--start-group                              \
 *              $@{MKL_LIB_PATH@}/libmkl_gf_lp64.a        \
 *              $@{MKL_LIB_PATH@}/libmkl_gnu_thread.a     \
 *              $@{MKL_LIB_PATH@}/libmkl_core.a           \
 *         -Wl,--end-group                                \
 *         -lgomp -lpthread"
 * ./configure --with-blas="$MKL" --with-lapack="$MKL"
 *  
 */
#ifdef USE_MKL_FIX
#include <stdlib.h>
#include <stdint.h>
void *_gfortran_internal_malloc64(uint64_t sz)
{
    return(calloc(1,(size_t)sz));
}
void *_gfortran_internal_malloc(uint32_t sz)
{
    return(calloc(1,(size_t)sz));
}
void _gfortran_internal_free(void *p)
{
    free(p);
}
#endif
