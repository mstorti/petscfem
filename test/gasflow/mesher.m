## usage: [xnod,icone,mesh]= mesher (XNOD,ICONE,H)
function [xnod,icone,mesh]= mesher (XNOD,ICONE,H)

  maxnode = rows(XNOD);
  
  NELEM = rows(ICONE);
  edge2elem = [(1:NELEM)';
  	       (1:NELEM)';
  	       (1:NELEM)';
  	       (1:NELEM)'];
  edges = [ICONE(:,[1 2]);
	   ICONE(:,[2 3]);
	   ICONE(:,[3 4]);
	   ICONE(:,[4 1])];
  edges = (min(edges')+maxnode*(max(edges')-1))';
  all_edges = edges;
  edge2elem = [edges edge2elem];

  flip = find(H(:,1)>H(:,2));
  if size(H,2)==3
    H(flip,:) = H(flip,[2 1 3]);
  else
    H(flip,:) = H(flip,[2 1 3 6 5 4]);
  endif

  edges = H(:,1:2);
  H(:,1:2) = [];
  edges = (min(edges')+maxnode*(max(edges')-1))';
  
  nedges = maxnode^2;
  HH = zeros(nedges,size(H,2));
  HH(edges,:) = H;

  set_flag=zeros(nedges,1);
  aux = create_set(edges);
  set_flag(aux) = ones(size(aux));
  clear aux

  front = edges;
  while length(front)>0
    edge = front(1);
    front(1)=[];
    ## printf("takes front %d %d\n",mesher_edge2(edge,maxnode));
    elems = find(edge2elem(:,1)==edge);
    elems = edge2elem(elems,2);
    for elem=elems'
      [HH,set_flagn] = \
      mesher_propagate_h(HH,set_flag, \
			 ICONE(elem,[1 2]),ICONE(elem,[4 3]));
      [HH,set_flagn] = \
      mesher_propagate_h(HH,set_flagn, \
			 ICONE(elem,[1 4]),ICONE(elem,[2 3]));
      front = [front; find(set_flagn& !set_flag)];
      set_flag = set_flagn;
    endfor
  endwhile

  xnod = []; icone=[];
  mesh.XNOD = XNOD;
  mesh.ICONE = ICONE;
  mesh.HH = HH;

  mesh.block_start = zeros(rows(icone));
  for k=1:rows(ICONE)
    [xx,ii] = mesher_isomap(mesh);
    nnod = rows(xnod);
    mesh.block_start(k) = nnod;
    xnod = [xnod; xx];
    icone = [icone; ii+nnod];
  endfor

  mesh.maxnode = maxnode;
  mesh.xnod = xnod;
  mesh.icone = icone;
  mesh.edge2elem = edge2elem;

  nnod = rows(xnod);
  nodemap = (1:nnod)';
  for edge=all_edges'
    elems = find(edge2elem(:,1)==edge);
    elems = edge2elem(elems,2);
    if length(elems)==1; continue; endif
    length(elems)==2 || error("edge connected to more than two elements");
    nodes1 = mesher_bound(mesh,edge,elems(1));
    nodes2 = mesher_bound(mesh,edge,elems(2));
    merr(xnod(nodes1,:)-xnod(nodes2,:)) < 1e-10 || \
	error("Not mactching edges");
    nodes = [nodes1 nodes2]';
    from = max(nodes)';
    to = min(nodes)';
    nodemap(from) = to;
  endfor

  while 1
    nodemap2 = nodemap(nodemap);
    elim = sum(nodemap!=nodemap2);
    if elim==0; break; endif
    ## printf("eliminated %d\n",elim);
    nodemap = nodemap2;
  endwhile
  
  ## printf("mapped %d\n",sum(nodemap!=(1:nnod)'));
  nnod2 = sum(nodemap==(1:nnod)');
  mask = (nodemap==(1:nnod)');
  map2 = cumsum(mask).*mask;

  nodemap = map2(nodemap);
  for k=1:columns(icone)
    icone(:,k) = nodemap(icone(:,k));
  endfor
  xnod = xnod(find(mask),:);

  mesh.maxnode = maxnode;
  mesh.xnod = xnod;
  mesh.icone = icone;
  mesh.XNOD = XNOD;
  mesh.ICONE = ICONE;
  mesh.edge2elem = edge2elem;
  mesh.HH = HH;
  mesh.nodemap = nodemap;

endfunction
