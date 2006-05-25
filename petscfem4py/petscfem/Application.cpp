// $Id: Application.cpp,v 1.1.2.2 2006/05/25 00:30:57 dalcinl Exp $

#include "Application.h"


PYPF_NAMESPACE_BEGIN

Application::~Application() 
{ 
  PYPF_DECREF(this->domain);
}

Application::Application()
  : Object(),
    domain(NULL)
{ }

Application::Application(const Application& app)
  : Object(app),
    domain(app.domain)
{
  PYPF_INCREF(this->domain);
}

Application::Application(Domain& domain)
  : Object(domain.getComm()),
    domain(&domain)
{
  PYPF_INCREF(this->domain);
}

Domain&
Application::getDomain() const
{
  return *this->domain;
}

PYPF_NAMESPACE_END
