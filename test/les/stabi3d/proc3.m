## proc2
nu = 0.000133;

start = 200;
last = columns(u);
uav = sum(u(1:Ny+1,start:last)')'/(last-start+1);

dudy = uav(2)/y(2);
tauw = nu*dudy;
ustar = sqrt(tauw);

U_bulk = sum(xcent(uav).*diff(y));

Cf = tauw/(0.5*U_bulk^2);
