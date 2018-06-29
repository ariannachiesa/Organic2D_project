close all

V = Poisson.V;

figure(1)
plot(V,'-b','LineWidth',3)
grid on
xlabel('Mesh Nodes')
ylabel('Potential [V]')

n = Poisson.n;

figure(2)
plot(n,'-r','LineWidth',3)
grid on
xlabel('Mesh Nodes')
ylabel('Electron density [m^-^3]')

resnrm = Poisson.res;

figure(3)
plot(resnrm,'-r','LineWidth',3)
grid on
title('Residual norm')

niter = Poisson.niter;

figure(4)
plot(niter,niter,'ob','MarkerFaceColor','b')
grid on
title('N. Iterations')