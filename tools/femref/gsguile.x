/* cpp arguments: -Wno-deprecated -I/home/mstorti/PETSC/petsc-2.1.6/include -I/home/mstorti/PETSC/petsc-2.1.6/bmake/linux -I/usr/local/mpich-1.2.5.2/include -I/usr/include/glib-1.2 -I/usr/lib/glib/include -I/home/mstorti/SOFT/NEWMAT/src -I/home/mstorti/SOFT/metis-4.0/Lib -I../.. -I/home/mstorti/SOFT/meschach-1.2 -DUSE_ANN -I/usr/local/ann_0.2/include -DUSE_SUPERLU -I/home/mstorti/SOFT/SuperLU -DUSE_DLEF -DUSE_SSL -I/home/mstorti/SOFT/SSL -DUSE_DX -DUSE_PTHREADS gsguile.cpp */
 scm_c_define_gsubr (s_comp_matrices_w, 6, 1, 0, (SCM (*)()) comp_matrices_w);
 scm_c_define_gsubr (s_fem_smooth_w, 8, 0, 0, (SCM (*)()) fem_smooth_w);
 scm_c_define_gsubr (s_elem2nod_proj_w, 4, 0, 0, (SCM (*)()) elem2nod_proj_w);
