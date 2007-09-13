// $Id$

#ifndef PF4PY_DOMAIN_H
#define PF4PY_DOMAIN_H

#include <memory>
#include <utility>
#include <vector>
#include "petscfem4py.h"

#include "Comm.h"
#include "Object.h"
#include "Options.h"
#include "DTable.h"
#include "Elemset.h"
#include "Mesh.h"
#include "Amplitude.h"
#include "Dofset.h"
#include "AppCtx.h"

PF4PY_NAMESPACE_BEGIN

class Domain 
  : public Object
{
  friend class Mesh;
  friend class Dofset;

private:
  Domain();
  Domain(const Domain& domain);
  Domain& operator=(const Domain& domain);

protected:
  MPI_Comm       comm;
  
  int            ndim;
  int            nnod;
  int            ndof;

  std::string    type;
  RefVal<AppCtx> appctx;
  RefVal<Mesh>   mesh;
  RefVal<Dofset> dofset;
  
  std::pair<VecScatter,Vec> scatter;
  
public:
  ~Domain();
  Domain(const std::string& type, int ndim, int nnod, int ndof);
  Domain(const std::string& type, int ndim, int nnod, int ndof, MPI_Comm comm);

  Comm getComm() const { return this->comm; }
  int  getNDim() const { return this->ndim; }
  int  getNNod() const { return this->nnod; }
  int  getNDof() const { return this->ndof; }

  const std::string& getType() const;

  Options& getOptions() const;
  void     setOptions(const Options& options);

  DTable<double>& getNodedata() const;
  void            setNodedata(const DTable<double>& nodedata);

  DTable<double>& getField(const std::string& name) const;
  void            setField(const std::string& name,
			   DTable<double>& data);

  Elemset& getElemset(int index) const;
  void     addElemset(const Elemset& elemset);
  
  void setFixation(int nn, const int    nodes[],
		   int nf, const int    fields[],
		   int nv, const double values[]);
  void setFixation(int nn, const int    nodes[],
		   int nf, const int    fields[],
		   Amplitude& amplitude);
  void setFixation(int nn, const int    nodes[],
		   int nf, const int    fields[],
		   int nv, const double values[],
		   Amplitude& amplitude);

  void setPeriodic(int n1, const int nodes1[],
		   int n2, const int nodes2[]);
  void setPeriodic(int n1, const int nodes1[],
		   int n2, const int nodes2[],
		   int nf, const int fields[]);
  void setConstraint(int nc, const double coeffs[],
		     int nn, const int    nodes[],
		     int nf, const int    fields[]);
  
  Mesh&   getMesh()   const;
  Dofset& getDofset() const;

  void setUp();

public:
  void allocateSolution(Vec& u) const;
  void allocateState   (Vec& x, const std::string& vec_type="") const;
  void allocateResidual(Vec& r, const std::string& vec_type="") const;
  void allocateJacobian(Mat& J, const std::string& mat_type="") const;

public:
  void assemble(double t, Vec x, 
		Vec r, Mat J) const;
  void assemble(double t1, Vec x1, double t0, Vec x0,
		Vec r, Mat J, double alpha=1.0) const;

  void assemble(const std::string& jobname,
		double t, Vec x,
		Vec r, Mat J) const;
  void assemble(const std::string& jobname,
		double t1, Vec x1, double t0, Vec x0,
		Vec r, Mat A, double alpha=1.0) const;

public:
  int  buildState   (double t, Vec u, Vec x) const;
  void buildSolution(double t, Vec x, Vec u) const;
  
};

PF4PY_NAMESPACE_END

#endif // PF4PY_DOMAIN_H

// Local Variables:
// mode: C++
// End:
