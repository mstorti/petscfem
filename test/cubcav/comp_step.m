source("./data.m.tmp");
if 1
  ## steps_dir = "/cdrom/STEPS-2003-SEP-01";
  steps_dir = "/cdrom/STEPS-2003-SEP-15b";
  pf_start_step = 116;
  refine = 1;
else
  steps_dir = "./STEPS";
  pf_start_step = 0;
  refine = 1;
endif
case_name = "cubcav";

dx_step = sscanf(getenv("dx_step"),"%d");
dx_step = dx_step 
rest = rem(dx_step,refine);
sss = pf_start_step + floor(dx_step/refine);
alpha = rest/refine;

tic;
u1 = aload([steps_dir "/" case_name ".state_" int2str(sss) ".tmp"]);
if alpha>0
  u2 = aload([steps_dir "/" case_name ".state_" int2str(sss+1) ".tmp"]);
  u = (1-alpha) * u1 + alpha * u2;
  printf("Using %f*u(%d)+%f*u(%d)\n",(1-alpha),sss,alpha,sss+1);
else
  u = u1;
  printf("Using u(%d)\n",sss);
endif
printf("Loading states %dsecs\n",toc);

tic;
u = pf_smooth(u,1,5);
printf("Smoothing %dsecs\n",toc);

asave(["./" case_name ".dx-state.tmp"],u);
