petsc_data_name = input("Enter Octave data file: > ","s");
system(["petscload.pl < ",petsc_data_name]);
source("tmp_petsc_data_script.m");
unlink("tmp_petsc_data_script.m");
clear petsc_data_name
