T=linspace(0,40,100);

figure(1);
plot(T,sw_dens(0*ones(size(T)),T,0),'-k');
set(gca,'fontsize',20);
xlabel('Temperature ($^\circ$C)','fontsize',20,'interpreter','latex');
ylabel('Density (kg/m$^3$)','fontsize',20,'interpreter','latex');

hold all;
plot(T,sw_dens(30*ones(size(T)),T,0),'-b');
grid on
grid minor

legend('Fresh water','Sea water (30PSU)')