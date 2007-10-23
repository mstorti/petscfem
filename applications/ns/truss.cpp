//__INSERT_LICENSE__
//$Id mstorti-v6-branch-1.0.2-9-g5c6b966 Tue Oct 23 16:52:48 2007 -0300$

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/texthf.h>
#include "./nsi_tet.h"
#include "./adaptor.h"
#include "./truss.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void truss::init() {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void truss::element_connector(const FastMat2 &xloc,
                              const FastMat2 &state_old,
                              const FastMat2 &state_new,
                              FastMat2 &res,FastMat2 &mat) {

}
