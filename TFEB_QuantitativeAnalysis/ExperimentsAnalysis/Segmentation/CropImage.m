%%                              CROPPING
% Sara Napolitano
% 8/3/2017

% -------------------------------------------------------------------------
% � una funzione che fa il crop della prima immagine acquisita al
% microscopio e restituisce il riquadro del crop
% -------------------------------------------------------------------------

function CropRect=CropImage(file)
    Im=imread(file);    % Im � un array contenente i dati dell'immagine
    i=imadjust(Im); % questa funzione aumenta il contrasto dell'immagine i
    [Icrop,CropRect]=imcrop(i);   % funzione che realizza il crop
                                  % Icrop � l'immagine dopo il cropping
                                  % CropRect � il riquadro del crop
    clear Im i Icrop;
end