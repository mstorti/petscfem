u1=aload("save.state.metis.tmp");
#u2=aload("save.state.hitchhiking.tmp");
#u3=aload("save.state.nearest_neighbor.tmp");
u4=aload("save.state.random.tmp");

tol = 1e-8;

#err = merr(u2-u1);
#printf("Hitchhiking partitioning OK ? > %d, [error: %g]\n",err<tol,err);

#err = merr(u3-u1);
#printf("Nearest neighbor partitioning OK ? > %d, [error: %g]\n",err<tol,err);

err = merr(u4-u1);
printf("Random partitioning OK ? > %d, [error: %g]\n",err<tol,err);

