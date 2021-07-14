%% Main_Fitting.m
%%% JULY, 2021

% WARNING: Please, note that an optimization problem can give little
% differences in results depending on the machine on which you run it, on
% the initial conditions you give etc.

warning off
clear all, close all, clc

%% Experimental Data

load('Data_TorinVariablePulses');

out = (out-min(out))./(max(out)-min(out));
expdata1 = iddata(out',in,Ts*60);

clear out in Ts

load('Data_TorinPulses');

out = (out-min(out))./(max(out)-min(out));
expdata2 = iddata(out',in,Ts*60);

clear out in Ts

expdata = merge(expdata1,expdata2);
expdata.exp = {'Variable Pulses';'Classical Pulses'};

clear expdata1 expdata2

%% Build NL model

% parameters: initial condition for the optimization problem
k = 1/70;
beta = 1/216;
a1 = beta;
b1 = beta;
c1 = 10.*k;
a3 = 1;
b3 = beta;
a = .00015;
b = .000088;

Parameters = {k,beta,a1,b1,c1,a3,b3,a,b}; % parameters to estimate
Order = [1,1,7];
InitialStates = [1 0 0 1 0 0 0]';
Ts = 0;

sys = idnlgrey('mTORmodel_feedback',Order,Parameters,InitialStates,Ts,'TimeUnit','s');
sys.SimulationOptions.Solver='ode45';

%% Estimation

opt = nlgreyestOptions('Display', 'on');
opt.SearchMethod='auto'; 
% opt.SearchOptions.Tolerance=1e-7; 
opt.SearchOptions.MaxIterations = 1000;

estModel = nlgreyest(expdata, sys, opt);    % estimated model

save('FittedData');

%% Figures for the Thesis
Ts = 15;

for k = 1:2
    
    esp = getexp(expdata,k);
    
    [simul,fit] = compare(esp,estModel);
    y = simul.OutputData;
    fit = round(fit.*100)./100;
    
    in = esp.InputData;
    out = esp.OutputData;
    dimEXP = length(in);
    ExpName = esp.exp{1};
    
    F = figure('Position', [1 1 720 435], 'DefaultAxesFontSize', 12, ...
                        'DefaultAxesLineWidth', 2.5, 'Renderer', 'Painters');

        subplot(3,1,1:2),hold on,
            plot([0:dimEXP-1].*Ts, out, 'LineWidth', 2.5, 'Color', [0.18,0.46,0.15]);
            plot([0:dimEXP-1].*Ts, y, 'LineWidth', 2.50, 'Color', [.30,0.75,0.93]);
            ylabel('Nuclear TFEB (%)')
            title('TFEB nuclear translocation')
            legend('Exp. Data',strcat('mTOR model (= ',num2str(fit),'%)'))
            set(gca, 'XLim', [0,length(in)*15-15], 'YLim', [0,1], 'XTick',...
                0:100:length(in)*15, 'XTickLabel', {},'Box','off', 'YTick',...
                0:.2:1, 'YTickLabel', 0:20:100);
            
        subplot(3,1,3)
            plot([0:dimEXP-1].*Ts, in, 'LineWidth', 2.5, 'Color', [.99 .55 .38]);
            ylabel('Input')
            xlabel('Time (min)'); 
            set(gca, 'XLim', [0,length(in)*Ts-Ts], 'XTick',  0:100:length(in)*Ts,... 
                'XTickLabel', {0:100:length(in)*Ts}, 'YLim', [-.1 1.1], 'Box', 'off',...
                'YTick', [0 1], 'YTickLabel', {'-Torin1','+Torin1'});

    print(F,strcat(ExpName),'-dsvg')
    print(F,strcat(ExpName),'-dpng')
    savefig(F,strcat(ExpName,'.fig'))

    close

end