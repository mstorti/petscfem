source("./data.m.tmp");
if 0
  ## steps_dir = "/cdrom/STEPS-2003-SEP-01";
  steps_dir = "./STEPS-2003-SEP-01";
  pf_start_step = 101;
else
  steps_dir = "./STEPS";
  pf_start_step = 0;
endif
refine = 5;
case_name = "cubcav";

cache.file_pattern = "STEPS/cubcav.state_%d.tmp";
cache.size = 5;
system("make fifo");
fifo = "comp_server.fifo";
fid = fopen(fifo,"r");
if fid==-1
  printf("couldn't open fifo");
  return;
endif

while 1

  printf("Waiting for step number at %s\n",fifo);
  [dx_step,count] = fscanf(fid,"%d");
  if count!=1; 
    printf("not found integer on fifo!!\n");
    continue;
  endif
  rest = rem(dx_step,refine);
  sss = pf_start_step + floor(dx_step/refine);
  alpha = rest/refine;

  [u1,cache] = comp_step_server(cache,sss);
  if alpha>0
    [u2,cache] = comp_step_server(cache,sss+1);
    u = (1-alpha) * u1 + alpha * u2;
    printf("Using %f*u(%d)+%f*u(%d)\n",(1-alpha),sss,alpha,sss+1);
  else
    u = u1;
    printf("Using u(%d)\n",sss);
  endif

  u = pf_smooth(u,1,5);

  asave(["./" case_name ".dx-state.tmp"],u);

endwhile
