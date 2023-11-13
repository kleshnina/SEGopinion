% Fisheries (get the equilibrium by letting system converge)

%close all

tic
c1 = .3; 
c2 = .2;
c11 = c1;
c12 = c2;
c21 = c1;
c22 = c2;
b = 12;
sigmam = 1.5;
sigmam2 = sigmam^2;
sigmaov = 1/sigmam2;
delta = .7; 
sigma2 = 1;
p1 = .3*10^(-5);
p2 = .9*10^(-5);
sigmanew = 1;
sigmau = .7;


m1v = 0:0.0025:2.2;
m2v = 0:0.0025:2.2; 

%m1v = 0:0.0005:2.2;
%m2v = 0:0.0005:2.2; 

Ustar = zeros(length(m1v), length(m2v));
X2 = Ustar;
X1 = Ustar;
U1 = Ustar;
U2 = Ustar;
Imb = Ustar;
Imbrel = Ustar;
Q = Ustar;

for i=1:length(m1v)
    m1 = m1v(i);
    if mod(i,50)==0
        disp(strcat('done with ', num2str((i-1)/881)));
    end
    for j=1:length(m2v)
        m2 = m2v(j);

        %%%%%%%%%%%% Get u*(m)=u*(m1, m2)
        
        x1t = [60];
        x2t = [20];
        u1t = [6];
        u2t = [2];

        for t=1:5000 %500000
            ft = b/(1+u1t);
            ht = m1*exp(-(u1t-m2)^2/sigmam^2)*(1-delta*exp(-u2t^2/sigmanew^2));
            st = u2t/(1+u2t);
            psit = 1/(1+c21*x1t+c22*x2t);
            sigmat = 1/(1+c11*x1t+c12*x2t);
            x1new = ft*psit*x2t;
            x2new = (1-ht)*st*sigmat*x1t;
            R0t = ft*(1-ht)*st*psit*sigmat;
            DR01t = -(b*(1 - exp(-((-m2 + u1t)^2/sigmam^2))*(1 - delta*exp(-(u2t^2/sigma2^2)))*m1)*u2t)/((1 +u1t)^2*(1 + u2t))*sigmat*psit+...
                     (2*b*exp(-((-m2 + u1t)^2/sigmam^2))*(1 - delta*exp(-(u2t^2/sigma2^2)))*m1*u2t*(u1t-m2))/((1 +u1t)*(1 + u2t)*sigmam^2)*sigmat*psit;
            
            u1new = u1t+sigmau*1/2*DR01t/R0t;
            
            DR02t = -(b*(1 - exp(-((-m2 + u1t)^2/sigmam^2))*(1 - delta*exp(-(u2t^2/sigma2^2)))*m1)*u2t)/((1 +u1t)*(1 + u2t)^2)*sigmat*psit+...
                    +b*(1 - exp(-((-m2 + u1t)^2/sigmam^2))*(1 - delta*exp(-(u2t^2/sigma2^2)))*m1)/((1 +u1t)*(1 + u2t))*sigmat*psit-...
                     (2*b*delta*exp(-((-m2 + u1t)^2/sigmam^2))*exp(-(u2t^2/sigma2^2))*m1*u2t^2)/((1 +u1t)*(1 + u2t)*sigma2^2)*sigmat*psit;
            u2new = u2t+sigmau*1/2*DR02t/R0t;
            
            x1t = x1new;
            x2t = x2new;
            u1t = u1new;
            u2t = u2new;
        
        end
        X2(i,j) = min(max(0.00001,x2t), 10000); %optimization across realistic values
        X1(i,j) = min(max(0.00001,x1t), 10000); %optimization across realistic values
        U1(i,j) = max(0,u1t);
        U2(i,j) = max(0,u2t);
        Imb(i,j) = min(X2(i,j)/X1(i,j),100);
        Imbrel(i,j) = min(abs(X1(i,j)-X2(i,j))/X1(i,j), 100);


        ft = b/(1+u1t);
        ht = m1*exp(-(u1t-m2)^2/sigmam^2)*(1-delta*exp(-u2t^2/sigmanew^2));
        st = u2t/(1+u2t);
        psit = 1/(1+c21*x1t+c22*x2t);
        sigmat = 1/(1+c11*x1t+c12*x2t);
        Q(i,j) = max(0,ht*st*psit*X1(i,j)-p1*m1-p2*m2);
        
    end
end


fig1 = figure;
imagesc(flipud(Q));
colorbar

xticklabels = 0:0.4:2.2;
xticks = linspace(1, size(X1, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

yticklabels = sort(0:0.4:2.2, 'desc');
yticks = linspace(1, size(X1, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)

ylabel('meshsize');
xlabel('effort');
title('Fisher''s profit');

box on
set(gca,'FontSize',14)
saveas(fig1, 'Fig1+2traits.png');
saveas(fig1, 'Fig1+2traits.fig');




fig11 = figure;
Qoff = mean(mean(Q));
Q(Q> 3*Qoff) = 3*Qoff;
imagesc(flipud(Q));
colorbar

xticklabels = 0:0.4:2.2;
xticks = linspace(1, size(X1, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

yticklabels = sort(0:0.4:2.2, 'desc');
yticks = linspace(1, size(X1, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)

ylabel('meshsize');
xlabel('effort');
title('Fisher''s profit');

box on
set(gca,'FontSize',14)
saveas(fig11, 'Fig11+2traits.png');
saveas(fig11, 'Fig11+2traits.fig');



fig2 = figure;
% % smooth out the plot
Imboff = mean(mean(Imb));
Imb(Imb> 3*Imboff) = 3*Imboff;
% 
imagesc(flipud(Imb));
colorbar

xticklabels = 0:0.4:2.2;
xticks = linspace(1, size(X1, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

yticklabels = sort(0:0.4:2.2, 'desc');
yticks = linspace(1, size(X1, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)

ylabel('meshsize');
xlabel('effort');
title('Imbalance index');

box on
set(gca,'FontSize',14)
saveas(fig2, 'Fig2+2traits.png');
saveas(fig2, 'Fig2+2traits.fig');




fig3=figure;

Imbreloff = mean(mean(Imbrel));
Imbrel(Imbrel> 3*Imbreloff) = 3*Imbreloff;

imagesc(flipud(Imbrel));
colorbar

xticklabels = 0:0.4:2.2;
xticks = linspace(1, size(X1, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

yticklabels = sort(0:0.4:2.2, 'desc');
yticks = linspace(1, size(X1, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)

ylabel('meshsize');
xlabel('effort');
title('Imbalance index');

box on
set(gca,'FontSize',14)

saveas(fig3, 'Fig3+2traits.png');
saveas(fig3, 'Fig3+2traits.fig');

toc

