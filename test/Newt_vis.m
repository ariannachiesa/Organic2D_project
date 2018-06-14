%% Script for Newton results visualization

close all

V1 = Newt.V;
n1 = Newt.n;

indexV = 1:2:(2*2403);
indexn = 2:2:(2*2403);
indexF = (2*2403+1):(2*2403+4);
indexI = (2*2403+5):(2*2403+6);

Vres = Res.res(indexV);
nres = Res.res(indexn);
Fres = Res.res(indexF);
Ires = Res.res(indexI);

figure(1)
plot(Vres,'-b','LineWidth',3)

figure(2)
plot(nres,'-r','LineWidth',3)

figure(3)
plot(Fres,'-ok','LineWidth',3)

figure(4)
plot(Ires,'-og','LineWidth',3)

Vout = V1 + Vres;
nout = n1 + nres;

figure(5)
plot(Vout,'-b','LineWidth',3)
title('Output potential of DD problem');
xlabel('[m]')
ylabel('[m]')
zlabel('potential [V]')

figure(6)
plot(nout,'-r','LineWidth',3)
title('Output electron density of DD problem');
xlabel('[m]')
ylabel('[m]')
zlabel('n [m^-^2]')