%% mTORmodel_feedback

function [xDOT,y] = mTORmodel_feedback(t,x,in,k,beta,a1,b1,c1,a3,b3,a,b,varagin)
    
    u = x(2)./(x(1)+x(2))-x(7);
    
    k1 = k.*(1-u);
    k_1 = k.*u;
    
    xDOT(1,1) = a1 - b1.*x(1) - c1.*x(3).*x(1); % mTOR_A
    xDOT(2,1) = c1.*x(3).*x(1) - b1.*x(2);  % mTOR_I
    xDOT(3,1) = a3.*in - b3.*x(3) + b1.*x(2) - c1.*x(3).*x(1);  % T_I
    
    xDOT(4,1) = k1.*(1-x(4)-x(5)-x(6))- k_1.*x(4) + beta.*x(6);
    xDOT(5,1) = k_1.*x(6) - k1.*x(5) + beta.*(1-x(4)-x(5)-x(6));
    xDOT(6,1) = k1.*x(5) - k_1.*x(6) - beta.*x(6);
    
    xDOT(7,1) = -a.*x(7)+ b.*(x(5)+ x(6));
    
    y = x(5) + x(6);
    
end