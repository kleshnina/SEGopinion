% Stackelberg equilibrium
% two-drug model

function odes = fdyn_pest(t,x,m1,m2)

% xvec = [x(1) x(2)];

%%%%%%%%% Parameter values
r = 0.7; K = 10000;
k1 = 1; k2 = 3; b1  = 5; b2 = 10; 
s1 = 0.01; s2 = 0.01;
% m1 = 0.9; 
% m2 = 1-m1;

%%%%%% Defining the system of ODEs

ode1 = x(1)*(r*((1-x(2))*(1-x(3))*K-x(1))/K-m1/(k1+b1*x(2))-m2/(k2+b2*x(3)));
ode2 = s1*(-r*(1-x(3))+(m1*b1)/(k1+b1*x(2))^2);
ode3 = s2*(-r*(1-x(2))+(m2*b2)/(k2+b2*x(3))^2);

odes = [ode1; ode2; ode3];

end
