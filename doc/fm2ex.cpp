  Chrono chrono;
  FastMat2 x(2,3,2),a(1,2),b(1,2),J(2,2,2);
  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);
  // Compute area of elements
  chrono.start();
  for (int ie=0; ie<nelem; ie++) {
    FastMat2::reset_cache();
    for (int k=1; k<=3; k++) {
      int node = ICONE(ie,k-1);
      x.ir(1,k).set(&XNOD(node-1,0)).rs();
    }
    x.rs();
    a.set(x.ir(1,2));
    a.rest(x.ir(1,1));

    b.set(x.ir(1,3));
    b.rest(x.ir(1,1));

    J.ir(1,1).set(a);
    J.ir(1,2).set(b);
    
    double area = J.rs().det()/2.;
    total_area += area;
    if (ie==0) {
      minarea = area;
      maxarea = area;
    }

    if (area>maxarea) maxarea=area;
    if (area<minarea) minarea=area;
  }
  printf("total_area %g, min area %g,max area %g, ratio: %g\n",
	 total_area,minarea,maxarea,maxarea/minarea);
  printf("Total area OK? : %s\n",
	 (fabs(total_area-1)<1e-8 ? "YES" : "NOT"));
  double cpu = chrono.elapsed();
  FastMat2::print_count_statistics();
  printf("CPU: %g, number of elements: %d\n"
	 "rate: %g [sec/Me], %g Mflops\n",
	 cpu,nelem,cpu*1e6/nelem,
	 nelem*FastMat2::operation_count()/cpu/1e6);
  FastMat2::deactivate_cache();
  FastMat2::void_cache();
