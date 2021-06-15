%%                              CROPPING
% Sara Napolitano
% 8/3/2017

% -------------------------------------------------------------------------
% È una funzione che fa il crop della prima immagine acquisita al
% microscopio e restituisce il riquadro del crop
% -------------------------------------------------------------------------

function CropRect=CropImage(file)
    Im=imread(file);    % Im è un array contenente i dati dell'immagine
    i=imadjust(Im); % questa funzione aumenta il contrasto dell'immagine i
    [Icrop,CropRect]=imcrop(i);   % funzione che realizza il crop
                                  % Icrop è l'immagine dopo il cropping
                                  % CropRect è il riquadro del crop
    clear Im i Icrop;
end