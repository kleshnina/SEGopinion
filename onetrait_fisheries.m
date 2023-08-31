% Cushing model_me (only one trait but let's consider two different stategies)

%moreruns image assumed p=20, costmesh=.4, costeffort=.2, and 10k runs
%until convergence 
syms x1 x2 u 

c11 = 0.01;
c12 = 0.005;
c21 = 0.001; 
c22 = 0.005;
b = 8;
sigmam = 0.3;
%effort = .4; 
%meshsize = .3;




%%%% Profit function
price = 20;
costmesh = 0.4;
costeffort = 0.2; 



sigma = 1/(1+c11*x1+c12*x2);
varphi = 1/(1+c21*x1+c22*x2);




meshv = [];
effortv = [];
Qv = [];
x1v = [];
x2v = [];
uv = [];
meshv1 = 0:0.01:1;
effv1 = 0:0.01:1;
Qm = zeros(length(meshv1), length(effv1));
Divx = Qm;
for mesht = 1:length(meshv1)

    for effortt = 1:length(effv1)
        meshsize = meshv1(mesht);
        effort = effv1(effortt);
        
        s = u/(1+u); %natural survival
        h = effort*exp(-10*(u-meshsize)^2);
        hs = (1-h);    %harvest survival (1-harvest)
        beta = 1/(1+u);
        f = b*beta;
        R0 = f*s*hs*varphi*sigma;
        DR0 = diff(R0, u);
        ratio = DR0/R0;
        x1s = f*varphi*x2;
        x2s = s*hs*sigma*x1;


        unew = 0.8;
        x1new = 50;
        x2new = 18;
        if subs(R0, [x1, x2, u], [x1new, x2new, unew])==0
            ratnew = 0;
        else
            ratnew = double(subs(ratio, [x1, x2, u], [x1new, x2new, unew]));
        end

        for j=1:10000
          x1new = double(subs(x1s, [x1, x2, u], [x1new, x2new, unew]));
          x2new = double(subs(x2s, [x1, x2, u], [x1new, x2new, unew]));
  
          unew = unew+sigmam*1/2*ratnew;
          
        if subs(R0, [x1, x2, u], [x1new, x2new, unew])==0
            ratnew = 0;
        else
            ratnew = double(subs(ratio, [x1, x2, u], [x1new, x2new, unew]));
        end

        end 
        %assume that after 500 iterations, the process converged
        hnew = double(subs(h, u, unew));
        snew = double(subs(s, u, unew));
        sigmanew = double(subs(sigma, [x1, x2], [x1new, x2new]));


        Qv = [Qv, hnew*snew*sigmanew*x1new*price-costmesh*meshsize-costeffort*effort];
        x1v = [x1v, x1new];
        x2v = [x2v, x2new];
        uv = [uv, unew];
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


% csvwrite('Qm_Cushing_me2_detail', Qm);
% csvwrite('Divm_Cushing_me2_detail', Divx);
% csvwrite('X1_Cushing_me2_detail', x1v);
% csvwrite('X2_Cushing_me2_detail', x2v);
% csvwrite('U_Cushing_me2_detail', uv);



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
