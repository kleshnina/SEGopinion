% Cushing model_me with two traits u1 and u2 

syms x1 x2 u1 u2 

c11 = 0.01;
c12 = 0.005;
c21 = 0.001; 
c22 = 0.005;
b = 8;
sigmam = 0.3;
%effort = .4; 
%meshsize = .3;




%%%% Profit function
price = 15;
costmesh = 1.2;
costeffort = .2; 



sigma = 1/(1+c11*x1+c12*x2);
varphi = 1/(1+c21*x1+c22*x2);




meshv = [];
effortv = [];
Qv = [];
x1v = [];
x2v = [];
u1v = [];
u2v = [];
meshv1 = 0:0.1:1;
effv1 = 0:0.1:1;
Qm = zeros(length(meshv1), length(effv1));
Divx = Qm;
for mesht = 1:length(meshv1)

    for effortt = 1:length(effv1)
        meshsize = meshv1(mesht);
        effort = effv1(effortt);
        
        s = u2/(1+u2); %natural survival
        h = effort*exp(-10*(u1-meshsize)^2)*(1-.5*exp(-u2^2));
        hs = (1-h);    %harvest survival (1-harvest)
        beta = 1/(1+u1);
        f = b*beta;
        R0 = f*s*hs*varphi*sigma;
        DR01 = diff(R0, u1);
        ratio1 = DR01/R0;
        DR02 = diff(R0,u2);
        ratio2 = DR02/R0;
        x1s = f*varphi*x2;
        x2s = s*hs*sigma*x1;


        u1new = 0.8;
        u2new = 0.7;
        x1new = 50;
        x2new = 18;
        if subs(R0, [x1, x2, u1, u2], [x1new, x2new, u1new, u2new])==0
            rat1new = 0;
            rat2new = 0;
        else
            rat1new = double(subs(ratio1, [x1, x2, u1, u2], [x1new, x2new, u1new, u2new]));
            rat2new = double(subs(ratio2, [x1, x2, u1, u2], [x1new, x2new, u1new, u2new]));
        end
        
        for j=1:10000
          x1new = double(subs(x1s, [x1, x2, u1, u2], [x1new, x2new, u1new, u2new]));
          x2new = double(subs(x2s, [x1, x2, u1, u2], [x1new, x2new, u1new, u2new]));
  
          u1new = u1new+sigmam*1/2*rat1new;
          u2new = u2new+sigmam*1/2*rat2new;
          
        if subs(R0, [x1, x2, u1, u2], [x1new, x2new, u1new, u2new])==0
            rat1new = 0;
            rat2new = 0;
        else
            rat1new = double(subs(ratio1, [x1, x2, u1, u2], [x1new, x2new, u1new, u2new]));
            rat2new = double(subs(ratio2, [x1, x2, u1, u2], [x1new, x2new, u1new, u2new]));

        end

        end 
        %assume that after 500 iterations, the process converged
        hnew = double(subs(h, [u1, u2], [u1new, u2new]));
        snew = double(subs(s, [u1, u2], [u1new, u2new]));
        sigmanew = double(subs(sigma, [x1, x2], [x1new, x2new]));


        Qv = [Qv, hnew*snew*sigmanew*x1new*price-costmesh*meshsize-costeffort*effort];
        x1v = [x1v, x1new];
        x2v = [x2v, x2new];
        u1v = [u1v, u1new];
        u2v = [u2v, u2new];
        meshv = [meshv, meshsize];
        effortv = [effortv, effort];
        
        Qm(mesht, effortt) = Qv(end);
        Divx(mesht, effortt) = (x1new-x2new)/x1new;   %relative difference between x1 and x2
        

    end



end



% stackelberg equilibrium: 
% let's see where there is a (local) maximum 
% locmaxidx = [];
% for kk = 2:length(Qv)-1
%     if Qv(kk)>Qv(kk-1) && Qv(kk)>Qv(kk+1)
%         locmaxidx = [locmaxidx, kk];
%     end
% 
% end
% 
% if length(locmaxidx)>0
%     % there is at least one local max
%     disp('the profit is');
%     disp(Qv(locmaxidx));
%     disp('the Stackelberg meshsize is');
%     disp(meshv(locmaxidx));
%     disp('the Stackelberg effort is');
%     disp(effortv(locmaxidx));
%     disp('the Stackelberg equilibrium x1 is');
%     disp(x1v(locmaxidx));
%     disp('the Stackelberg equilibrium x2 is');
%     disp(x2v(locmaxidx));
%     disp('the Stackelberg equilibrium u is');
%     disp(uv(locmaxidx));
% 
% end


% meshvector = 0:.1:1;
% Qm = zeros(length(meshvector), length(meshvector));
% 
% for i1=1:length(Qv)
%     effi = mod(i1, length(meshvector))+1;
%     meshi = floor(i1/length(meshvector))+1;
%     Qm(meshi, effi) = Qv(i1);
% end

%csvwrite('Qm2', Qm);

Qmflip = flipud(Qm);

figure(1)
imagesc(Qmflip);
colorbar
xticks([1:1:11]);
xticklabels(meshv1);
xlabel('effort');
yticks(1:1:11);
%yticklabels(0.8:0.02:1);
yticklabels(sort(effv1, 'desc'));
ylabel('meshsize');
title('Q');


Divflip = flipud(Divx);
Divflip(isnan(Divflip))=1;
Divflip(Divflip<0)=0;

figure(2)
imagesc(Divflip);
colorbar
xticks([1:1:11]);
xticklabels(meshv1);
xlabel('effort');
yticks(1:1:11);
%yticklabels(0.8:0.02:1);
yticklabels(sort(effv1, 'desc'));
ylabel('meshsize');
title('Imbalance Index');
