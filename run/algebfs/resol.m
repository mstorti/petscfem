###key resol.m

#[u, errorA, total_iters] = gmres(0*b, b, "fullsys",[1e-10 200 0]);
[p, errorH, total_iters] = gmres(0*bH, bH, "redsys",[1e-10 200 0]);
#[p, errorH, total_iters] = bicgstab(0*bH, bH, "redsys",[1e-10 200 0]);
