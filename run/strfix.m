#usage: s = strfix (a)
function s = strfix (a)

  s = setstr(a+32*(a==0));

endfunction
