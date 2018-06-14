%% Script for Newton results visualization

close all

figure(1)
plot(V,'-b','LineWidth',3)
hold on
plot(Vold,'--b','LineWidth',3)

figure(2)
plot(n,'-r','LineWidth',3)
hold on
plot(nold,'--r','LineWidth',3)

figure(3)
plot(F,'-k','LineWidth',3)
hold on
plot(Fold,'--ok','LineWidth',3)

figure(4)
plot(I,'-g','LineWidth',3)
hold on
plot(I,'--og','LineWidth',3)