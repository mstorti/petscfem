// $Id$

#ifndef PF4PY_APPCTX_H
#define PF4PY_APPCTX_H

#include <memory>
#include "petscfem4py.h"
#include "Object.h"

PF4PY_NAMESPACE_BEGIN

class AppCtx
  : public Object
{
 
  friend class Domain;

protected:
  AppCtx& operator=(const AppCtx&);

#if !defined(SWIG)
public:
  struct Args {
  public:
    typedef ::arg_list Impl;
  protected:
    Args(const Args&);
    Args& operator=(const Args&);
  protected:
    std::auto_ptr< Args::Impl > impl;
  public:
    virtual ~Args() = 0;
    Args();
    virtual Impl* argl() const { return this->impl.get(); }
    virtual const char* job()  const = 0;
    virtual const Time* time() const = 0;
  };
protected:
  std::auto_ptr<Args> args;
#endif

protected:
  virtual void setup(Domain*) { }
protected:
  static  void assemble(const Domain* domain, const Args* args);
protected:
  virtual bool sdgraph(const Domain* domain,
		       std::vector<int>& dofs,
		       std::vector<int>& xadj,
		       std::vector<int>& adjncy) const { return false; }
  virtual bool profile(const Domain* domain,
		       std::vector<int>& xadj,
		       std::vector<int>& adjncy) const { return false; }

  virtual void assemble(const Domain* domain,
			double t, Vec x,
			Vec r, Mat J) const 
  { this->assemble(domain,"",t,x,r,J); }
  virtual void assemble(const Domain* domain,
			double t1, Vec x1,
			double t0, Vec x0, 
			Vec r, Mat J, double alpha) const 
  { this->assemble(domain,"",t1,x1,t0,x0,r,J,alpha); }

  virtual void assemble(const Domain* domain,
			const std::string& jobname,
			double t, Vec x,
			Vec r, Mat J) const = 0;
  virtual void assemble(const Domain* domain,
			const std::string& jobname,
			double t1, Vec x1,
			double t0, Vec x0, 
			Vec r, Mat J, double alpha) const = 0;
public:
  ~AppCtx(); 
  AppCtx(const AppCtx& app);
  AppCtx();
};

PF4PY_NAMESPACE_END



PF4PY_NAMESPACE_BEGIN

class AppNS
  : public AppCtx
{
private:
  AppNS& operator=(const AppNS&);

protected:
  virtual bool sdgraph(const Domain* domain,
		       std::vector<int>& dofs,
		       std::vector<int>& xadj,
		       std::vector<int>& adjncy) const;
  virtual bool profile(const Domain* domain,
		       std::vector<int>& xadj,
		       std::vector<int>& adjncy) const;

  virtual void assemble(const Domain* domain,
			const std::string& jobname,
			double t, Vec x,
			Vec r, Mat J) const;
  virtual void assemble(const Domain* domain,
			const std::string& jobname,
			double t1, Vec x1,
			double t0, Vec x0, 
			Vec r, Mat J, double alpha) const;
public:
  ~AppNS();
  AppNS(const AppNS& ns);
  AppNS();

};

PF4PY_NAMESPACE_END


#endif // PF4PY_APPCTX_H

// Local Variables:
// mode: C++
// End:
