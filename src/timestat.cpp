//__INSERT_LICENSE__
//$Id: timestat.cpp,v 1.1 2001/11/01 19:11:33 mstorti Exp $

#include <cmath>
#include <cstdio>

#include <src/timestat.h>

void TimeStat::init() {
  histo = new int[nbin+1];
  cumul = new double[nbin+1];
  for (int j=0; j<nbin+1; j++) {
    histo[j]=0;
    cumul[j]=0;
  }
  out = &histo[nbin];
  cout = &cumul[nbin];
  dt = (log(t_max)-log(t_min))/nbin;
}

void TimeStat::add(double t) {
  if (t<t_min || t >= t_max) {
    *out++;
    *cout += t;
  } else {
    int jbin = int(floor((log(t)-log(t_min))/dt));
    histo[jbin]++;
    cumul[jbin] += t;
  }
}

void TimeStat::print_stat() {
  int total = 0;
  double ctotal=0.;
  for (int j=0; j<nbin+1; j++) {
    total += histo[j];
    ctotal += cumul[j];
  }
  if (*out>0) 
    printf("warning: entries outside range %d (%g \%) , %g secs (%g\%)\n",
	   (*out), (*out)/total, (*cout) / ctotal);
  for (int j=0; j<nbin; j++) {
    printf("%g <= t < %g : %d (%g\%), %g secs (%g\%)\n",
	   t_min*exp(j*dt),
	   t_min*exp((j+1)*dt),
	   histo[j],
	   double(histo[j])/double(total),
	   cumul[j],
	   cumul[j]/ctotal
	   );
  }
}

// TimeStat time_stat;
