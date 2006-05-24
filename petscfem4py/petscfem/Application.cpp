// $Id: Application.cpp,v 1.1.2.1 2006/05/24 21:02:33 dalcinl Exp $

#include "Application.h"


PYPF_NAMESPACE_BEGIN

Application::~Application() 
{ 
  PYPF_DECREF(this->domain);
}

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
