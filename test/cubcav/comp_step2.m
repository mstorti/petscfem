source("./data.m.tmp");
if 1
  ## steps_dir = "/cdrom/STEPS-2003-SEP-01";
  steps_dir = "./STEPS-2003-SEP-01";
  pf_start_step = 101;
else
  steps_dir = "./STEPS";
  pf_start_step = 0;
endif
refine = 5;
case_name = "cubcav";

cache.file_pattern = [steps_dir "/cubcav.state_%d.tmp"];
cache.size = 5;
system("make fifo");
fifo = "comp_server.fifo";
fifo2 = "comp_server.fifo2";

while 1

  fid = fopen(fifo,"r");
  if fid==-1
    printf("couldn't open fifo");
    return;
  endif
  printf("Waiting for step number at %s\n",fifo);
  [dx_step,count] = fscanf(fid,"%d");
  fclose(fid);
  if count!=1; 
    printf("not found integer on fifo!!\n");
    continue;
  endif
  printf("Received step %d\n",dx_step);
  if dx_step<0; break; endif
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

  u = pf_smooth(u,1,4);

  asave(["./" case_name ".dx-state.tmp"],u);

  
  fid = fopen(fifo2,"w");
  if fid==-1
    printf("couldn't open fifo2");
    return;
  endif
  fprintf(fid,"OK\n");
  fclose(fid);

endwhile

system("make clean_fifo");
