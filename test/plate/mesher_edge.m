#usage: 
function [edge2,invert] = mesher_edge(mesh,edge)

  maxnode = mesh.maxnode;
  edge2 = (min(edge')+maxnode*(max(edge')-1))';
  invert = edge(:,1)>edge(:,2);

endfunction
