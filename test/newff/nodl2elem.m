%usage: function  elems = nodl2elem (nodl)
%      returns a list of elemes from a node list (1d mesh)
function  elems = nodl2elem (nodl)

  nnod=rows(nodl);
  elems=[nodl(1:nnod-1) nodl(2:nnod)];

endfunction
