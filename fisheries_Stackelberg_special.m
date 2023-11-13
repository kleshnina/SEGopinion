% Fisheries  (assuming c11=c21 and c12=c21 ->psi = sigma
close all

c1 = .3; 
c2 = .2; 
b = 12;
sigmam = 1.5;
sigmam2 = sigmam^2;
sigmaov = 1/sigmam2;
p1 = .3*10^(-5);
p2 = .9*10^(-5);

%m1v = 0:0.02:1.2;
%m2v = 0:0.02:2.2;
m1v = 0:0.00005:1.2;
m2v = 0:0.00005:2.2; 

Ustar = zeros(length(m1v), length(m2v));
X2 = Ustar;
X1 = Ustar;
Q = Ustar;
Imb = Ustar;
Imbrel = Imb;

for i=1:length(m1v)
    m1 = m1v(i);
    for j=1:length(m2v)
        m2 = m2v(j);

        %%%%%%%%%%%% Get u*(m)=u*(m1, m2)
        ustarmf =@(u) exp((m2-u)^2*sigmaov)*sigmam2*(u-1)-m1*sigmam2*(u-1)+2*m1*(m2-u)*u*(1+u);
        
        ustarm = fzero(ustarmf, rand*10);
 %        ustarm = fzero(ustarmf, [0, 50]);

        count = 0;
        while ustarm <=0 && count<100
            count = count+1;
            ustarm = fzero(ustarmf, rand*10*count);
        end
        if ustarm<=0 
            disp('ustarm not found');
            Ustar(i,j) = 1000*ustarm;

        else
            Ustar(i,j) = ustarm;
        end
        
        %%%%%%%%%%% Get x*(m,u*(m)) 
        fu = b/(1+ustarm);
        hu = m1*exp(-(ustarm-m2)^2/sigmam^2);
        su = ustarm/(1+ustarm);
        gu = fu*(-1 + hu)*su;

        %%% get x1, x2 analytically
        if hu<1
            x1eq = (sqrt(fu*(1-hu)*su)-1)/(c1+c2*sqrt((1-hu)*su/fu));
            x2eq = sqrt((1-hu)*su/fu)*x1eq;
            
            X2(i,j) = max(0,x2eq);
            X1(i,j) = max(0,x1eq);
            Imb(i,j) = X2(i,j)/X1(i,j);
            Imbrel(i,j) = abs(X1(i,j)-X2(i,j))/X1(i,j);
        end


        %%%%%%%%%%% Get Q matrix (profit) income = h*s*psi*x1
        Q(i,j) = max(0,hu/(1-hu)*X2(i,j)-p1*m1-p2*m2);
    end
end

figure(1)
imagesc(flipud(Q));
colorbar

xticklabels = 0:0.2:1.2;
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


figure(2)
imagesc(flipud(Imb));
colorbar

xticklabels = 0:0.2:1.2;
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




figure(3)
imagesc(flipud(Imbrel));
colorbar

xticklabels = 0:0.2:1.2;
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


