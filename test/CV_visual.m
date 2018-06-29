%% CV curve output visualization

close all

C = CVcurve.C;
V = CVcurve.V;

figure(1)
plot(V,C,'-b','LineWidth',3)
xlabel('Gate Voltage [V]');
ylabel('Capacitance [F]');
grid on