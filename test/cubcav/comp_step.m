source("./data.m.tmp");
## steps_dir = "./STEPS-2003-SEP-01";
if 0
  steps_dir = "/cdrom/STEPS-2003-SEP-01";
  pf_start_step = 101;
else
  steps_dir = "./STEPS";
  pf_start_step = 0;
endif
refine = 5;
case_name = "cubcav";

dx_step = sscanf(getenv("dx_step"),"%d");
dx_step = dx_step 
rest = rem(dx_step,refine);
sss = pf_start_step + floor(dx_step/refine);
alpha = rest/refine;

u1 = aload([steps_dir "/" case_name ".state_" int2str(sss) ".tmp"]);
if alpha>0
  u2 = aload([steps_dir "/" case_name ".state_" int2str(sss+1) ".tmp"]);
  u = (1-alpha) * u1 + alpha * u2;
  printf("Using %f*u(%d)+%f*u(%d)\n",(1-alpha),sss,alpha,sss+1);
else
  u = u1;
  printf("Using u(%d)\n",sss);
endif

u = pf_smooth(u,1,5);

asave(["./" case_name ".dx-state.tmp"],u);
