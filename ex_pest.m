tspan=[0,250];
x0 = [1000; 0.5; 0.5];
Y = 0.8; a = 1; c1 = 0.3; c2 = 0.3; K=10000;
m=0.01:0.1:0.99;
pop=zeros(length(m),length(m));
u1=zeros(length(m),length(m));
u2=zeros(length(m),length(m));
Q=zeros(length(m),length(m));

for i=1:length(m)
    for j=1:length(m)
        m1=m(i);
        m2=m(j);
        [t,x]=ode45(@(t,x) fdyn_pest(t,x,m1,m2),tspan,x0);
        l=length(t);
        pop(i,j) = x(l,1);
        u1(i,j) = x(l,2);
        u2(i,j) = x(l,3);
        Q(i,j) = Y*(1-a*(pop(i,j)/K)^2)-c1*m1-c2*m2;
    end
end

%%%%%%%% Plotting equilibria 

h = round(length(m)/2);
l = length(m);

subplot(2,2,1)

poppl = flipud(pop);

image(poppl,'CDataMapping','scaled')
colormap parula
caxis([0 10000])
colorbar

title('Population size')
xlabel('Pesticide 1')
ylabel('Pesticide 2')
set(gca,'fontsize',14)

yticks([1 h l])
yticklabels({'1','0.5','0'})
xticks([1 h l])
xticklabels({'0','0.5','1'})

subplot(2,2,3)

u1pl = flipud(u1);

image(u1pl,'CDataMapping','scaled')
colormap parula
u1min = round(min(min(u1)),1);
u1max = round(max(max(u1)),1);
caxis([u1min u1max])
colorbar

title('Resistance rate against Pesticide 1')
xlabel('Pesticide 1')
ylabel('Pesticide 2')
set(gca,'fontsize',14)

yticks([1 h l])
yticklabels({'1','0.5','0'})
xticks([1 h l])
xticklabels({'0','0.5','1'})

subplot(2,2,4)

u2pl = flipud(u2);

image(u2pl,'CDataMapping','scaled')
colormap parula
u2min = round(min(min(u2)),1);
u2max = round(max(max(u2)),1);
caxis([u2min u2max])
colorbar

title('Resistance rate against Pesticide 2')
xlabel('Pesticide 1')
ylabel('Pesticide 2')
set(gca,'fontsize',14)

yticks([1 h l])
yticklabels({'1','0.5','0'})
xticks([1 h l])
xticklabels({'0','0.5','1'})

subplot(2,2,2)

Qpl = flipud(Q);

image(Qpl,'CDataMapping','scaled')
colormap parula
Qmin = round(min(min(Q)),1);
Qmax = round(max(max(Q)),1);
caxis([Qmin Qmax])
colorbar

title('Farmers profit')
xlabel('Pesticide 1')
ylabel('Pesticide 2')
set(gca,'fontsize',14)

yticks([1 h l])
yticklabels({'1','0.5','0'})
xticks([1 h l])
xticklabels({'0','0.5','1'})

%%%%%%%% Plotting population dynamics
% 
% subplot(2,1,1)
% 
% plot(t,x(:,1),'Color',[0.4660 0.6740 0.1880],'LineWidth',4) %green
% 
% box off
% % title('Population dynamics')
% xlabel('Time')
% ylabel('Population size')
% set(gca,'fontsize',16)
% % legend('sensitive','resistant')
% ylim([0 10000])
% yticks([0 5000 10000])
% % xticks([0 50 100 150 200 250])
% 
% subplot(2,1,2)
% 
% plot(t,x(:,2),'Color',[0.4940 0.1840 0.5560],'LineWidth',4) %purple
% hold on
% plot(t,x(:,3),'Color',[0.3010 0.7450 0.9330],'LineWidth',4) %blue
% 
% box off
% % title('Resistance rate evolution')
% xlabel('Time')
% ylabel('Resistance rate')
% set(gca,'fontsize',16)
% % legend('drug 1','drug 2')
% ylim([-0.1 1])
% yticks([0 0.5 1])