%%                          VERIFICA IMMAGINI
% Sara Napolitano
% 10/3/2017

% -------------------------------------------------------------------------
% � una funzione che mi restituisce 1 se trova l'immagine in memoria
% -------------------------------------------------------------------------

function s=ExistImage(file)
    if exist(file)==2   % quando � uguale a 2 esiste un file con il nome che viene passato a exist
        s=1;    % se esiste mi restituisce 1
    else
        s=0;    % se non esiste mi restituisce 0
    end
end