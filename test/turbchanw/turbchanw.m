rho=1000;			# kg/m^3
g=9.8;				# m/sec^2
nu=1;				# m^2/sec
y_wall_plus=8.;			# adim.
von_Karman_cnst=1.;		# 
roughness=0.;			# .01 m
A0=-0.217;			# adim
A1=5.50;			# adim
A2=2.50;			# adim
A3=8.5;				# adim
Chezy=110;			# ????




h=0.1;				# m
slope = 1/10000;

## Water viscosity, from (Heat, mass and Momentum Transfer, W. Rosenhow
## & H. Choi. Prentice-Hall 1961
##    T[F]      nu[ft^2/hr]       T[C]      nu[m^2/sec]
##     32           0.070         0.00       1.75e-06
##     60           0.044        15.55       1.10e-06
##    100           0.027        37.77       6.75e-07
##    200           0.012        93.33       3.00e-07

## Mean flow

q=sqrt(Chezy^2*hin*slope)
Fr=Ust/sqrt(g*hst)

Pk=g*q^3/C^2;
epsilon=Pk/h;

Pe=C2*sqrt(Cmu)*g^(5/4)*q^4/(h*sqrt(D)*Chezy^2/5);
k=C2*epsilon^2*h/Pe;

