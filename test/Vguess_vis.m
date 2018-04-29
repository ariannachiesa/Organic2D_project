% Vguess visualization

tx = unique(msh.p(1,:));
ty = unique(msh.p(2,:));

tz = msh.f(:)';

%n = ceil(length(tz)/length(xx));
n = length(tx);

rows = length(tx);
cols = length(tz)/length(tx);

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
plot(xx,yy);
