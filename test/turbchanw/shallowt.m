##__INSERT_LICENSE__
## $Id: shallowt.m,v 1.1 2003/02/17 14:53:21 mstorti Exp $
sigma_k=1.;			# Turbulence model constants
sigma_e=1.3;
C_mu = 0.09;
C_1=1.44;
C_2=1.92;
D = 1.;

Chezy = 110;			# Friction with the bottom

rho=1000;			# kg/m^3
gravity=9.8;			# m/sec^2
nu=1e-6;			# m^2/sec
y_wall_plus=8.;			# adim.
von_Karman_cnst=1.;		# 
roughness=0.;			# .01 m
A0=-0.217;			# adim
A1=5.50;			# adim
A2=2.50;			# adim
A3=8.5;				# adim
