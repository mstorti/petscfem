#!/usr/bin/perl -pi.bak

s/sles\.h/petscsles.h/g;
s/vec\.h/petscvec.h/g;
s/sles\.h/petscsles.h/g;
s/mat\.h/petscmat.h/g;
s/(\W)Scalar(\W)/$1PetscScalar$2/g;
s/(\W)Viewer(\W)/$1PetscViewer$2/g;
s/VIEWER/PETSC_VIEWER/g;
s/OptionsGetString/PetscOptionsGetString/g;
s/PETSC_VIEWER_FORMAT_ASCII_MATLAB/PETSC_VIEWER_ASCII_MATLAB/g;
