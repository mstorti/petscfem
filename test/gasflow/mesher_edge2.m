#usage: 
function edge2 = mesher_edge (mesh,edge)

  maxnode = mesh.maxnode;
  n2 = rem(edge,maxnode);
  n1 = round((edge-n2)/maxnode)+1;
  edge2 = [n1,n2];

endfunction
