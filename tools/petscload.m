##__INSERT_LICENSE__
## $Id: petscload.m,v 1.2 2003/01/08 15:54:26 mstorti Exp $
petsc_data_name = input("Enter Octave data file: > ","s");
system(["petscload.pl < ",petsc_data_name]);
source("tmp_petsc_data_script.m");
unlink("tmp_petsc_data_script.m");
clear petsc_data_name
