#usage: 
function edge2 = mesher_edge(mesh,edge)

  maxnode = mesh.maxnode;
  edge2 = (min(edge')+maxnode*(max(edge')-1))';

endfunction
