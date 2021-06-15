%% SEGMENTAZIONE OFF-LINE
% Sara Napolitano
% 21/06/19

% -------------------------------------------------------------------------
% Off-line segmentaition algoeithm for HeLa TFEB-GFP
% Run first FastER segmentation on phase contrast images
% -------------------------------------------------------------------------

close all, clear all, clc

colordef white

% =========================================================================
% Variabili di INPUT:
% =========================================================================

% Folder where the images are saved
Path = strcat(pwd,'ImagePath\');    % example path, it is empty

% experiment name
name_exp = 'HBSSPulses';

fcsuffix = 'c1.tif';
greensuffix = 'c2.tif';
redsuffix = 'c3.tif';
masksuffix = 'c1_m00_mask.png';

Ts = 15; 
fields = 1;     % how many chambers are acquired

dimEXP = 90;    % how many time points

% input
% Normal Pulses
in = [zeros(3*60/Ts,1);ones(3*60/Ts,1);zeros(3*60/Ts,1);ones(3*60/Ts,1);zeros(3*60/Ts,1);ones(dimEXP-5*3*60/Ts,1)];
% % No First Pulses
% in =[zeros(3*60/Ts,1);zeros(3*60/Ts,1);zeros(3*60/Ts,1);ones(3*60/Ts,1);zeros(3*60/Ts,1);ones(dimEXP-5*3*60/Ts,1)];
% % Single Pulse
% in = ones(dimEXP,1);
% % Variable Pulses 
% in = [zeros(60/Ts,1);ones(0.5*60/Ts,1);zeros(3*60/Ts,1);ones(60/Ts,1);zeros(3*60/Ts,1);ones(2*60/Ts,1);zeros(3*60/Ts,1);ones(3*60/Ts,1);zeros(3*60/Ts,1);ones(dimEXP-3*4*60/Ts-60/Ts,1)];
% % Single Pulses
% in = [zeros(60/Ts,1);ones(12*60/Ts,1);zeros(dimEXP-9*60/Ts,1)];

% -------------------------------------------------------------------------
% struct initialization
CELLS = cell(3,fields);

for i=1:fields
    CELLS{3,i} = 1;
end

% =========================================================================
% Cropping:
% =========================================================================
for i=1:fields
    clc
    tempstr = strcat('Cropping the immage from field number:',num2str(i));
    disp(tempstr);
    
    image = strcat('t',num2str(1,'%.2d'),'xy',num2str(i,'%.2d'),fcsuffix);
    CropRect = CropImage(strcat(Path,image));
    vn = genvarname(strcat('CropRect',num2str(i)));
    eval([vn,'=CropRect;']);
    
    save(strcat(vn,'.mat'),'CropRect');
    clear tempstr CropRect vn
end

clear i, close all, clc
save(strcat(name_exp,'.mat'));

% =========================================================================
% Time points cycle
% =========================================================================

for t = 1:dimEXP  
    
    disp(['ITR:' num2str(t) '/' num2str(dimEXP)])
    CELLS = trackingCELLS(CELLS,Path,fields,t);
    save(strcat(name_exp,'.mat'));
    
end