## usage: r = atv(x)
function r = latv (x)

  global linear_atv_data

  r = linear_atv_data.A * x;

endfunction
