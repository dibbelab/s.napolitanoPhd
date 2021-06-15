%% Feedback_Fun.m

function xDOT = Feedback_Fun(t,x)%,u)

global u
% -------------------------------------------------------------------------
% Parameters definition

k = 1./70; % de/phosphorylation rate assigned value
beta = 1./216; % transport rate assigned value
u_cap = u - x(4);
a = .00015;
b = .000088;

% -------------------------------------------------------------------------
% de/phosphorylation rates
c1 = k; 
c_1 = k;
c2 = k;
c_2 = k;
k1 = c1.*(1-u_cap);
k_1 = c_1.*(u_cap);
k2 = c2.*(1-u_cap);
k_2 = c_2.*(u_cap);
% transport rates
beta1 = beta;
beta2 = beta;

% -------------------------------------------------------------------------
% xDOT definition: 
% x(1) = [TFEBcp], x(2) = [TFEBn], x(3) = [TFEBnp]
% [TFEBc] = 1-x(1)-x(2)-x(3)

xDOT(1,1) = k1.*(1-x(1)-x(2)-x(3))- k_1.*x(1) + beta2.*x(3);
xDOT(2,1) = k_2.*x(3) - k2.*x(2) + beta1.*(1-x(1)-x(2)-x(3));
xDOT(3,1) = k2.*x(2) - k_2.*x(3) - beta2.*x(3);
xDOT(4,1) = -a.*x(4)+ b.*(x(2)+ x(3));