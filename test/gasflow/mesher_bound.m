#usage: 
function nodes = mesher_bound (mesh,edge,elem);

  if size(edge,2)>2 & size(edge,1)==1
    nodes = [];
    for k=1:size(edge,2)-1
      nodess = mesher_bound (mesh,edge([k k+1]));
      k==1 || nodes(length(nodes)) == nodess(1) || \
	  error("couldn't concatenate set od edges");
      if k!=1; nodess(1) = []; endif
      nodes = [nodes; nodess];
    endfor
    return;
  endif

  invert_flag = 0;
  if size(edge,2) == 2
    invert_flag = edge(1)>edge(2);
    edge = mesher_edge(mesh,edge);
  endif

  if nargin<=2
    elems = find(mesh.edge2elem(:,1)==edge);
    elems = mesh.edge2elem(elems,2);
    if !(length(elems)==1 || length(elems)==2)
      mesher_edge2(mesh,edge)
      error("edge should be connected to 1 or 2 elements");
    endif
    elem = elems(1);
  endif

  ICONE = mesh.ICONE(elem,:);

  edges = [ICONE([1 2]);
	   ICONE([2 3]);
	   ICONE([3 4]);
	   ICONE([4 1])];
  edges = mesher_edge(mesh,edges);
  k = find(edges==edge);
  m = mesh.HH(edges(1),1);
  n = mesh.HH(edges(2),1);
  if k==1
    nodes = 1+(0:m)'*(n+1);
    invert_flag = invert_flag != (ICONE(1)>ICONE(2));
  elseif k==2
    nodes = m*(n+1)+(1:n+1)';
    invert_flag = invert_flag != (ICONE(2)>ICONE(3));
  elseif k==3
    nodes = n+1+(0:m)'*(n+1);
    invert_flag = invert_flag != (ICONE(4)>ICONE(3));
  elseif k==4
    nodes = (1:n+1)';
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
