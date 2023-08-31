% Stackelberg equilibrium
% two-drug model

syms x0 x1 x2 m1 m2 u1 u2 

x = [x0 x1 x2];

%%%%%%%%% Parameter values
rmax = 0.45; d = 0.01; K = 10000;
k1 = 5; k2 = 5; b1  = 10; b2 = 10; g1 = 0.5; g2=0.5;
s1 = 1; s2 = 1;

a1 = 0.15; a2 = 0.9; a0 = 1; a3 = 0.9;
A = [a0 a1 a1; a2 a0 a3; a2 a3 a0];

Qmax = 1; c = 0.5; w1 = 0.5; w2 = 0.2; r1 = 0.4; r2 = 0.4;

%%%%%%% Fitness Functions
G0 = rmax*(1-(x*A(1,:)')/K)-d-m1/k1-m2/k2;
G1 = rmax*exp(-g1*u1)*(1-(x*A(2,:)')/K)-d-m1/(k1+b1*u1)-m2/k2;
G2 = rmax*exp(-g2*u2)*(1-(x*A(3,:)')/K)-d-m1/k1-m2/(k2+b2*u2);

%%%%%%% Derivatives of fitness Functions
dG1du1 = -g1*rmax*exp(-g1*u1)*(1-(x*A(2,:)')/K)+m1*b1/(k1+b1*u1)^2;
dG2du2 = -g2*rmax*exp(-g2*u2)*(1-(x*A(3,:)')/K)+m2*b2/(k2+b2*u2)^2;

%%%%%%% Leader function
Q = Qmax - c*((x0+x1+x2)/K)^2 - w1*m1^2 - w2*m2^2 - r1*u1^2 - r2*u2^2;

%%%%%%% Equilibria of Ecological Dynamics
x0eq = (K*((a0^2-a3^2)*(-d-m1/k1-m2/k2+rmax)+a1*(-a0+a3)*exp(g1*u1)*(-d-m2/k2+...
    exp(-g1*u1)*rmax-m1/(k1+b1*u1))+a1*(-a0+a3)*exp(g2*u2)*(-d-m1/k1+...
    exp(-g2*u2)*rmax-m2/(k2+b2*u2))))/((a0-a3)*(a0^2-2*a1*a2+a0*a3)*rmax);
x1eq = (K*(a2*(-a0+a3)*(-d-m1/k1-m2/k2+rmax)+(a0^2-a1*a2)*exp(g1*u1)*(-d-m2/k2+...
    exp(-g1*u1)*rmax-m1/(k1+b1*u1))+(a1*a2-a0*a3)*exp(g2*u2)*(-d-m1/k1+...
    exp(-g2*u2)*rmax-m2/(k2+b2*u2))))/((a0-a3)*(a0^2-2*a1*a2+a0*a3)*rmax);
x2eq = (K*(a2*(-a0+a3)*(-d-m1/k1-m2/k2+rmax)+(a1*a2-a0*a3)*exp(g1*u1)*(-d-m2/k2+...
    exp(-g1*u1)*rmax-m1/(k1+b1*u1))+(a0^2-a1*a2)*exp(g2*u2)*(-d-m1/k1+...
    exp(-g2*u2)*rmax-m2/(k2+b2*u2))))/((a0-a3)*(a0^2-2*a1*a2+a0*a3)*rmax);


%%%%%% Planar Map at Ecological Dynamics/Fitness Functions at Ecol. Equi. 
G0eq = subs(G0, [x0,x1,x2],[x0eq,x1eq,x2eq]);
G1eq = subs(G1, [x0,x1,x2],[x0eq,x1eq,x2eq]);
G2eq = subs(G2, [x0,x1,x2],[x0eq,x1eq,x2eq]);

dG1du1eq = subs(dG1du1, [x0,x1,x2],[x0eq,x1eq,x2eq]);
dG2du2eq = subs(dG2du2, [x0,x1,x2],[x0eq,x1eq,x2eq]);

%%%%%% Leader function at Ecol. Equi.
Qeq = subs(Q, [x0, x1, x2], [x0eq, x1eq, x2eq]);



% we need to find u1*,u2* (as a function of m1,m2) such that Gi is maximized
mleft = 0;
mright = 1;
mvec = mleft:0.1:mright;
msize = length(mvec);
U1starvec = zeros(msize,msize);
U2starvec = zeros(msize,msize);
for i = 1:msize
    for j = 1:msize
    %for each m value, we calculate u^*
    me1 = mvec(i);
    me2 = mvec(j);
    G1m = subs(dG1du1eq, [m1,m2], [me1, me2]);
    G2m = subs(dG2du2eq, [m1,m2], [me1, me2]);
    
    % we want (u1,u2) such that the partial G1m/partial u1 and partial
    % G2m/partial u2 is zero. Thus, we minimize abs(G1m_u1)+abs(G2m_u2)
    
    u = [u1, u2];
    f = matlabFunction(abs(G1m)+abs(G2m),'Vars',{u});
    
    %find for each u^* such that lamdbaim is maximized
    opts = optimset('fminsearch');
    opts.Display = 'iter';
    opts.TolX = 1.e-12;
    opts.MaxFunEvals = 1000;

    sol = fminsearchbnd(f, [0.1 0.1], [0 0], [1 1], opts);
    
    u1star = sol(1);
    u2star = sol(2);
    
    U1starvec(i,j) = u1star;
    U2starvec(i,j) = u2star;

    end
end

%mivec Uistarvec is the corresponding vectors where ui* gives the maximum
%values for corresponding mvec values.


%%%%%%%%%% STACKELBERG Equilibrium
% to calculate the Stackelberg equilibrium, we cycle through the array
% combinations [m1vec, m2vec, U1starvec, U2starvec] and check where the maximum Q
% occurs.
Qvec = zeros(msize,msize);
for i=1:msize
    for j=1:msize
    
    mmaxt1 = mvec(i);
    mmaxt2 = mvec(j);
    u1start = U1starvec(i,j);
    u2start = U2starvec(i,j);
    Qeqt = double(subs(Qeq, [m1, m2, u1, u2], [mmaxt1, mmaxt2, u1start, u2start]));
    
    Qvec(i,j) = Qeqt;
    end
end

[row,col] = find(Qvec == max(Qvec,[],"all"));
% idx = find(Qvec == max(Qvec));
if length(row)>1
    disp('there are more than one ecologically enlightened optimum');
end

% Stackelberg Equilibrium
disp('the Stackelberg equilibrium is at:');
disp(strcat('m_1= ', num2str(mvec(row))));
disp(strcat('m_2= ', num2str(mvec(col))));
disp(strcat('u_1= ', num2str(U1starvec(row,col))));
disp(strcat('u_2= ', num2str(U2starvec(row,col))));
disp(strcat('Qmax= ', num2str(Qvec(row,col))));

%%%% Assigning Stackelberg equilibrium values
EqS_m1 = mvec(row); EqS_m2 = mvec(col); EqS_u1 = U1starvec(row,col); EqS_u2 = U2starvec(row,col); EqS_Q = Qvec(row,col);
