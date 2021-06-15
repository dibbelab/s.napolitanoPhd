%% OpenLoop_Fun.m

function xDOT = OpenLoop_Fun(t,x)%,u)

global u
% -------------------------------------------------------------------------
% Parameters definition

k = 1./70; % de/phosphorylation rate assigned value
beta = 1./216; % transport rate assigned value

% -------------------------------------------------------------------------
% de/phosphorylation rates
c1 = k; 
c_1 = k;
c2 = k;
c_2 = k;
k1 = c1.*(1-u);
k_1 = c_1.*(u);
k2 = c2.*(1-u);
k_2 = c_2.*(u);
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