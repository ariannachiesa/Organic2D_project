% CV curve plot in mesh nodes with y = t_ins

close all

tx = unique(msh.p(1,:));

V = CVcurve.Vg;
C = CVcurve.C;

num = length(tx);
cont = 1;

for i = 1:num:length(C)
  c(cont) = C(i);
  cont++;
end

figure(1)
plot(V,c,'-b','LineWidth',2)