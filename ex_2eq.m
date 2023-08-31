tspan=[0,500];
x0 = [1000;1000; 0.1; 0.1];
% x0 = [1000;1000; 0.1];

[t,x]=ode45(@fdyn_2eq,tspan,x0);

xtot = x(:,1)+x(:,2);

subplot(2,1,1)

plot(t,x(:,1),'Color',[0.4660 0.6740 0.1880],'LineWidth',4) %green
hold on
plot(t,x(:,2),'Color',[0 0.4470 0.7410],'LineWidth',4)
% plot(t,xtot)

box off
% title('Population dynamics')
xlabel('Time')
ylabel('Population size')
set(gca,'fontsize',16)
% legend('sensitive','resistant')
ylim([0 10000])
yticks([0 5000 10000])
% xticks([0 50 100 150 200 250])

subplot(2,1,2)

plot(t,x(:,3),'Color',[0.4940 0.1840 0.5560],'LineWidth',4) %purple
hold on
plot(t,x(:,4),'Color',[0.3010 0.7450 0.9330],'LineWidth',4) %blue

box off
% title('Resistance rate evolution')
xlabel('Time')
ylabel('Resistance rate')
set(gca,'fontsize',16)
% legend('drug 1','drug 2')
ylim([-0.1 1])
yticks([0 0.5 1])