source("data.m.tmp");
ghia;

load -force sqcav.ny.tmp
u=aload("sqcav.lu.tmp");
usbp1=aload("sqcav.iisd_sbp1.tmp");
usbp2=aload("sqcav.iisd_sbp2.tmp");
usbp4=aload("sqcav.iisd_sbp4.tmp");
usbp8=aload("sqcav.iisd_sbp8.tmp");
usbp16=aload("sqcav.iisd_sbp16.tmp");

tol=1e-10;
erro = merr(usbp1-u);
printf("IISD iisd_subpart=1. error < tol OK ? %d (error = %g, tol = %g)\n", \
       erro<tol,erro,tol);

tol=1e-5;
erro = merr(usbp2-u);
printf("IISD iisd_subpart=2. error < tol OK ? %d (error = %g, tol = %g)\n", \
       erro<tol,erro,tol);

erro = merr(usbp4-u);
printf("IISD iisd_subpart=4. error < tol OK ? %d (error = %g, tol = %g)\n", \
       erro<tol,erro,tol);

erro = merr(usbp8-u);
printf("IISD iisd_subpart=8. error < tol OK ? %d (error = %g, tol = %g)\n", \
       erro<tol,erro,tol);

erro = merr(usbp16-u);
printf("IISD iisd_subpart=16. error < tol OK ? %d (error = %g, tol = %g)\n", \
       erro<tol,erro,tol);
