// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: mmove.h,v 1.6 2002/11/30 23:42:06 mstorti Exp $

#ifndef MMOVE_H
#define MMOVE_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class  mesh_move : public adaptor { 
private:
  FastMat2 G, J, dNdxi, xlocp, xloc0, res_Dir;
  //#define USE_NEWMAT
public: 
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
  double distor_fun(FastMat2 & xlocp);
  virtual void init_dfun() {}
  virtual double distor_fun_G(FastMat2 &G)=0;
};

class mesh_move_eig : public mesh_move {
private:
#ifdef USE_NEWMAT
  SymmetricMatrix  GG;
  DiagonalMatrix D;
#else
  FastMat2 D;
#endif
  double c_volume, c_distor, distor_exp;
public:
  void init_dfun();
  double distor_fun_G(FastMat2 & G);
};

class mesh_move_rcond : public mesh_move {
private:
  FastMat2 iG;
public:
  double distor_fun_G(FastMat2 & G);
};

class mmove_hook {
private:
  int write_mesh(const State &s,const char *filename,const int append=0);
  Dofmap *dofmap;
  Mesh *mesh;
public:
  void init(Mesh &mesh_a,Dofmap &dofmap_a,
	    TextHashTableFilter *options,const char *name) {
    mesh = &mesh_a; dofmap = &dofmap_a;
  }
  void time_step_pre(double time,int step) {}
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
};

#endif
