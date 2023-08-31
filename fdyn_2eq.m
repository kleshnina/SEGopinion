% Stackelberg equilibrium
% two-drug model

function odes = fdyn_2eq(t,x)

xvec = [x(1) x(2)];

%%%%%%%%% Parameter values
rmax = 0.45; d = 0.01; K = 10000;
k1 = 5; k2 = 5; b1  = 10; b2 = 10; g1 = 0.5; g2=0.5;
s1 = 1; s2 = 1;

% rmax = 0.45; d = 0.01; K = 10000;
% k1 = 4; k2 = 2; b1  = 10; b2 = 10; g1 = 0.8; g2=0.8;
% s1 = 1; s2 = 1;

%% new param
% rmax = 0.45; d = 0.01; K = 10000;
% k1 = 5; k2 = 5; b1  = 10; b2 = 10; g1 = 0.5; g2=0.5;
% s1 = 1; s2 = 1;

a1 = 0.15; a2 = 0.9;
A = [1 a1; a2 1];

m1 = 0.4; m2 = 0.5;

%%%%%% Defining the system of ODEs

ode1 = x(1)*(rmax*(1-(xvec*A(1,:)')/K)-d-m1/k1-m2/k2);
ode2 = x(2)*(rmax*exp(-g1*x(3)-g2*x(4))*(1-(xvec*A(2,:)')/K)-d-m1/(k1+b1*x(3))-m2/(k2+b2*x(4)));
ode3 = s1*(-g1*rmax*exp(-g1*x(3)-g2*x(4))*(1-(a2*x(1)+x(2))/K)+(m1*b1)/(k1+b1*x(3))^2);
ode4 = s2*(-g2*rmax*exp(-g1*x(3)-g2*x(4))*(1-(a2*x(1)+x(2))/K)+(m2*b2)/(k2+b2*x(4))^2);

odes = [ode1; ode2; ode3; ode4];

% ode1 = x(1)*(rmax*(1-(xvec*A(1,:)')/K)-d-m1/k1);
% ode2 = x(2)*(rmax*exp(-g1*x(3))*(1-(xvec*A(2,:)')/K)-d-m1/(k1+b1*x(3)));
% ode3 = s1*(-g1*rmax*exp(-g1*x(3))*(1-(a2*x(1)+x(2))/K)+(m1*b1)/(k1+b1*x(3))^2);
% 
% odes = [ode1; ode2; ode3];

end
