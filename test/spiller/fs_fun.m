#usage: 
function y = fs_fun(x)

  global spiller_data
  y = spline(spiller_data.xfs(:,1),spiller_data.xfs(:,2),x);

endfunction
