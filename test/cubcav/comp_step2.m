source("./data.m.tmp");
if 1
  steps_dir = "/cdrom/STEPS-2003-SEP-15b";
  pf_start_step = 116;
  refine = 10;
else
  steps_dir = "./STEPS";
  pf_start_step = 0;
endif
case_name = "cubcav";

cache.file_pattern = [steps_dir "/cubcav.state_%d.tmp"];
cache.size = 5;
system("make fifo");
fifo = "comp_server.fifo";
fifo2 = "comp_server.fifo2";

while 1

  start=time;
  printf("Waiting for step number at %s\n",fifo);
  tic;
  fid = fopen(fifo,"r");
  if fid==-1
    printf("couldn't open fifo");
    return;
  endif
  printf("Waiting %g\n",toc);
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

  if alpha>0
    printf("Using %f*u(%d)+%f*u(%d)\n",(1-alpha),sss,alpha,sss+1);
  else
    printf("Using u(%d)\n",sss);
  endif

  [u1,cache] = comp_step_server(cache,sss);
  if alpha>0
    [u2,cache] = comp_step_server(cache,sss+1);
    u = (1-alpha) * u1 + alpha * u2;
  else
    u = u1;
  endif
  printf("loaded states\n");

  tic;
  u = pf_smooth(u,0.25,4);
  printf("Smoothing %g\n",toc);

  tic;
  state_file = ["./" case_name ".dx-state.tmp"];
  if 0
    asave(state.tmp,u);
  else
    fid = fopen(state_file,"w");
    uu = u';
    [m,n] = size(uu);
    count = fwrite (fid,uu,"double");
    fclose(fid);
  endif;  
  printf("Saving %g\n",toc);
  
  fid = fopen(fifo2,"w");
  if fid==-1
    printf("couldn't open fifo2");
    return;
  endif
  fprintf(fid,"OK\n");
  fclose(fid);
  printf("Total %g\n",time-start);

endwhile

system("make clean_fifo");
