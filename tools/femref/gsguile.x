/* cpp arguments: -Wno-deprecated -I/u/mstorti/PETSC/petsc-2.1.6/include -I/u/mstorti/PETSC/petsc-2.1.6/bmake/linux -I/usr/local/mpich-1.2.5.2/include -I/usr/include/glib-1.2 -I/usr/lib/glib/include -I/u/mstorti/SOFT/NEWMAT/src -I/u/mstorti/SOFT/metis-4.0/Lib -I../.. -I/u/mstorti/SOFT/meschach-1.2 -DUSE_ANN -I/usr/local/ann_0.2/include -DUSE_SUPERLU -I/u/mstorti/SOFT/SuperLU -DUSE_DLEF -DUSE_SSL -I/u/mstorti/SOFT/SSL -DUSE_DX -DUSE_PTHREADS gsguile.cpp */
 scm_c_define_gsubr (s_comp_matrices_w, 6, 1, 0, (SCM (*)()) comp_matrices_w);
 scm_c_define_gsubr (s_elem2nod_proj_w, 6, 0, 0, (SCM (*)()) elem2nod_proj_w);
 scm_c_define_gsubr (s_nod2elem_proj_w, 4, 0, 0, (SCM (*)()) nod2elem_proj_w);
