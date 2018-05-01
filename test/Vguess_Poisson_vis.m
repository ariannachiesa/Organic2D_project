% Script for the visualization of the initial guess for Poisson (Vguess)
% and the outputs of Poisson

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

figure(2)
[xx,yy] = meshgrid(tx,ty);
plot(xx,yy);

figure(3)
mesh (tx, ty, vv);

figure(4)
mesh (tx, ty, nn);

figure(5)
plot(iter,res);
