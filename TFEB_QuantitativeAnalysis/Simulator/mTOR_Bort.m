%% Feedback_Fun.m

function xDOT = mTOR_Bort(t,x)%,u)

global u

% -------------------------------------------------------------------------
% Parameters definition
k = 0.0027; %1/70;
beta = 5.6e-04; %1/216;
a2 = 2e-04; %beta;
b2 = 0;  %0.97e-04;  %beta;
c1 = 2.16e-06; %10.*k;
a1 = 29.57; %1;
b1 = 0.0021; %beta;

a = 0.86e-04;  %.00015;
b = 0.89e-04;  %.000088;
    
u_cap = x(7)./(x(6)+x(7))-x(4);

% -------------------------------------------------------------------------
% de/phosphorylation rates
k1 = k.*(1-u_cap);
k_1 = k.*(u_cap);
% k2 = k.*(1-u_cap);
% k_2 = k.*(u_cap);

% transport rates
% beta1 = beta;
% beta2 = beta;

% -------------------------------------------------------------------------
% xDOT definition:

% TFEB
xDOT(1,1) = k1.*(1-x(1)-x(2)-x(3))- k_1.*x(1) + beta.*x(3);
xDOT(2,1) = k_1.*x(3) - k1.*x(2) + beta.*(1-x(1)-x(2)-x(3));
xDOT(3,1) = k1.*x(2) - k_1.*x(3) - beta.*x(3);

% Feedback
xDOT(4,1) = -a.*x(4)+ b.*(x(2)+ x(3));

% mTOR + Torin
xDOT(5,1) = a1.*u - b1.*x(5) + b2.*x(7) - c1.*x(5).*x(6);
xDOT(6,1) = a2 - b2.*x(6) - c1.*x(5).*x(6);
xDOT(7,1) = c1.*x(5).*x(6) - b2.*x(7);