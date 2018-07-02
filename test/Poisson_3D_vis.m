% Script for the visualization of the initial guess for Poisson (Vguess)
% and the outputs of Poisson

% N.B: load also file Vguess_visualization

close all

tx = unique(msh.p(1,:));
tz = msh.f(:)';

V = Poisson.V;
n = Poisson.n;
iter = Poisson.niter;
res = Poisson.res;

num = length(tx);

cont = 1;
for i = 1:num:length(msh.p(2,:))
  ty(cont) = msh.p(2,i);
  cont++;
end

rows = length(tz)/length(tx);
cols = length(tx);

zz = zeros(rows,cols);
vv = zeros(rows,cols);
nn = zeros(rows,cols);

for i = 1:num
  cont = 1;
  for j = i:num:length(tz)
    zz(cont,i) = tz(j);
    cont++;
  end
end

for i = 1:num
  cont = 1;
  for j = i:num:length(V)
    vv(cont,i) = V(j);
    nn(cont,i) = n(j);
    cont++;
  end
end

figure(1)
mesh (tx, ty, zz);
title('Initial Guess for Poisson problem');
xlabel('[m]')
ylabel('[m]')
zlabel('potential [V]')

figure(2)
[xx,yy] = meshgrid(tx,ty);
plot(xx,yy);
title('Scheme of the mesh');
xlabel('number of refinements along x-dimension')

figure(3)
mesh (tx, ty, vv);
title('Output potential of Poisson problem');
xlabel('[m]')
ylabel('[m]')
zlabel('potential [V]')

figure(4)
mesh (tx, ty, nn);
title('Output electron density of Poisson problem');
xlabel('[m]')
ylabel('[m]')
zlabel('n [m^-^2]')

if(iter ~= 1)
  for i = 1:iter
    niter(i) = i;
  end

  figure(5)
  plot(niter,res,'-or','LineWidth',3,'MarkerFaceColor','r');
  xlabel('Number of iterations')
  ylabel('Residual norm')
  
endif