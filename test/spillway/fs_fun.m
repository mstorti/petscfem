#usage: 
function y = fs_fun(x)

  global spillway_data
  y = spline(spillway_data.xfs(:,1),spillway_data.xfs(:,2),x);

endfunction
