function [xnod,icone,state]= proc1 (filename,ext)
%
%  [xnod,icone,state]= proc1 (filename,ext)
%

  eval(['load ' filename '.nod.' ext])
  eval(['xnod = ' filename ';'])
  eval(['load ' filename '.con.' ext])
  eval(['icone = ' filename ';'])
  eval(['load ' filename '.state.' ext])
  eval(['state = ' filename ';'])

end
