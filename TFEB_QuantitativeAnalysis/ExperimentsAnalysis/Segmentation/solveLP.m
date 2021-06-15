%% trackingYEA.m
%%% APRIL 14, 2017
% GIANSIMONE

function [M, STATUS] = solveLP(oldP, newP)

N1 = size(oldP,1);

N2 = size(newP,1);

COST = zeros(1,N1*N2);

Aeq = sparse(N2,N1*N2);


for i = 1:N1    % per tutti le cellule identificate nel ciclo prima
    
    for j = 1:N2    % per tutte le cellule identificate in questo ciclo
        
        COST(1,(((i-1)*N2)+j)) = norm(oldP(i,:) - newP(j,:));   % calcola la funzione di costo (distanza euleriana tra ogni vecchia cellula e ogni nuova cellula)
        
    end
    
end


for i = 1:N2    % per tutte le cellule identificate in questo ciclo
    
    Aeq(i,i:N2:N1*N2) = 1;  % cambia il vincolo (matrice diagonale di 1)

end


Beq = ones(1,N2);   % vincolo


LB = sparse(N1*N2,1);   % il vincolo inferiore è zero

UB = [];    % non c'è vincolo superiore

options = optimoptions('linprog','Algorithm','dual-simplex');   % by Sara
[xMIN, ~, STATUS] = linprog(COST, [], [], Aeq, Beq, LB, UB, options);    % funzione di minimizzazione del costo
                                                                % questa funzione di minimizzazione trova in quali casi le cellule corrispondono

M = reshape(xMIN, N2, N1)'; % cambio la dimensione alla matrice xMIN per renderla come mi serve

M = uint32(M);  % converto gli elementi dell'array in degli interi a 32 bit

end