#usage: 
function nodes = mesher_row (mesh,elem,edge,row);

  ICONE = mesh.ICONE(elem,:);

  edges = [ICONE([1 2]);
	   ICONE([2 3]);
	   ICONE([3 4]);
	   ICONE([4 1])];
  edges = mesher_edge(mesh,edges);
  if (size(edge,2)==2) 
    [edge,invert_flag] = mesher_edge(mesh,edge);
  endif

  k = find(edges==edge);
  m = mesh.HH(edges(1));
  n = mesh.HH(edges(2));
  if k==1
    nodes = row+(0:m)'*(n+1);
    invert_flag = invert_flag != (ICONE(1)>ICONE(2));
  elseif k==2
    nodes = (m-row+1)*(n+1)+(1:n+1)';
    invert_flag = invert_flag != (ICONE(2)>ICONE(3));
  elseif k==3
    nodes = n+2-row+(0:m)'*(n+1);
    invert_flag = invert_flag != (ICONE(4)>ICONE(3));
  elseif k==4
    nodes = (row-1)*(n+1)+(1:n+1)';
    invert_flag = invert_flag != (ICONE(1)>ICONE(4));
  endif

  if invert_flag
    nodes = nodes(length(nodes):-1:1);
  endif
  nodes = nodes + mesh.block_start(elem);

  if struct_contains(mesh,"nodemap")
    nodes = mesh.nodemap(nodes);
  endif

endfunction
