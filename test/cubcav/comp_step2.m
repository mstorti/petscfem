## cache.file_pattern = "STEPS-2003-SEP-01/cubcav.state_%d.tmp";
cache.file_pattern = "STEPS/cubcav.state_%d.tmp";
cache.size = 5;
system("mkfifo ");
cache.fifo_name = "comp_step_server.fifo";

while 1
  step = floor(rand*10)
  [u,cache] = comp_step_server(cache,step);
  pause
endwhile
