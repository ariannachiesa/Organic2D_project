% Vguess visualization

tx = unique(msh.p(1,:));
tz = msh.f(:)';

n = length(tx);

cont = 1;
for i = 1:n:length(msh.p(2,:))
  ty(cont) = msh.p(2,i);
  cont++;
end

rows = length(tz)/length(tx);
cols = length(tx);

zz = zeros(rows,cols);

for i = 1:n
  cont = 1;
  for j = i:n:length(tz)
    zz(cont,i) = tz(j);
    cont++;
  end
end

figure(1)
mesh (tx, ty, zz);

figure(2)
[xx,yy] = meshgrid(tx,ty);
plot(xx,yy);
