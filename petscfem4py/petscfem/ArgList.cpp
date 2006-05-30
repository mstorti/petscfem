// $Id: ArgList.cpp,v 1.1.2.1 2006/05/30 20:13:20 dalcinl Exp $

#include "ArgList.h"

#include <fem.h>


PYPF_NAMESPACE_BEGIN

ArgList::~ArgList() 
{  delete (ArgList::Base*)(*this); }

ArgList::ArgList()
  : Handle(new ArgList::Base)
{ }

PYPF_NAMESPACE_END
