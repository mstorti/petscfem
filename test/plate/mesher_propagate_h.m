#usage: 
function [HHn,set_flagn] = mesher_propagate_h (HH,set_flag,Edge1,Edge2);

  maxnode = sqrt(length(set_flag));
  flag1 = Edge1(1) > Edge1(2);
  flag2 = Edge2(1) > Edge2(2);
  edge1=min(Edge1)+maxnode*(max(Edge1)-1);
  edge2=min(Edge2)+maxnode*(max(Edge2)-1);
  
  if set_flag(edge1) & !set_flag(edge2)
    HH(edge2,:) = HH(edge1,:);
    set_flag(edge2)=1;
#      printf("Propagates refinement from %d %d to %d %d, refinement %d\n",
#  	   Edge1,Edge2,HH(edge1));
  elseif set_flag(edge2) & !set_flag(edge1)
    HH(edge1,:) = HH(edge2,:);
    set_flag(edge1)=1;
#      printf("Propagates refinement from %d %d to %d %d, refinement %d\n",
#  	   Edge2,Edge1,HH(edge1));
  endif

  HHn = HH;
  set_flagn = set_flag;

endfunction
