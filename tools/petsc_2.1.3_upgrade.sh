#!/usr/bin/perl -pi.bak

perl -pi.bak -e 's/sles\.h/petscsles.h/g'  \ 
 -e 's/vec\.h/petscvec.h/g;' \
 -e 's/sles\.h/petscsles.h/g;' \
 -e 's/mat\.h/petscmat.h/g;' \
 -e 's/(\W)Scalar(\W)/$1PetscScalar$2/g;' \
 -e 's/(\W)Viewer(\W)/$1PetscViewer$2/g;' \
 -e 's/VIEWER/PETSC_VIEWER/g;' \
 -e 's/OptionsGetString/PetscOptionsGetString/g;' \
 -e 's/PETSC_VIEWER_FORMAT_ASCII_MATLAB/PETSC_VIEWER_ASCII_MATLAB/g;' $*
