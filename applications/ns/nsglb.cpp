//__INSERT_LICENSE__
//$Id$
#include <src/debug.h>
#include <malloc.h>

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/sttfilter.h>
#include <src/pfmat.h>
#include <src/hook.h>
#include <src/iisdmatstat.h>

#include <applications/ns/nsi_tet.h>

vector<double>    data_pts;
vector<ElemToPtr> elemset_pointer;

int TSTEP           = 0;
int fractional_step = 0;
int reuse_mat       = 0;

double FLUID_TIME=0;

WallData wall_data;
