files = glob("outvector*.sal");
nfiles = rows(files);
ctimes = zeros(nfiles,1);

for k=1:nfiles
  info = stat(deblank(strfix(files(k,:))));
  ctimes(k) = info.ctime;
endfor
clear info

files = sortby(ctimes,files);

U = [];
for k = 1:nfiles
  UU = aload(deblank(strfix(files(k,:))));
  U = [U;UU];
endfor

clear UU

N = 861;
nstep = rows(U)/N;
if nstep != round(nstep)
  error("does not contain an integer number of vectors\n");
endif

hh = zeros(N,nstep);
for k=1:nstep
  hh(:,k) = U((k-1)*N+(1:N),4);
endfor

