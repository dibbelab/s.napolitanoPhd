%% Main_Simulator.m
%%% JUNE, 2021

warning off
clear all, close all, clc

% -------------------------------------------------------------------------
sims_list = {{'FIG4-6B', 'OpenLoop', [0.3 0.7], [0.7 0.3 0]};...
    {'FIG4-7B', 'Feedback', [0.5 1], [0.65 0.35 0 0.18]}};

dimEXP = 1290;

for z = 1:length(sims_list)
    
    ExpName = sims_list{z}{1};
    Mod = sims_list{z}{2};  % MOD := OpenLoop op Feedback
    InputLev = sims_list{z}{3};
    InitCond = sims_list{z}{4};
    
    % ---------------------------------------------------------------------
    % Simulation
    
    [Time,NucTFEB] = simulator(Mod,InputLev,InitCond,dimEXP);
    
    Ts = 15;
    in = [zeros(3*60/Ts,1);ones(3*60/Ts,1);zeros(3*60/Ts,1);ones(3*60/Ts,1);zeros(3*60/Ts,1);ones((dimEXP-5*3*60)/Ts,1)];
    % ---------------------------------------------------------------------
    % Plots
    F = figure('Position', [1 1 720 435], 'DefaultAxesFontSize', 12, ...
                'DefaultAxesLineWidth', 2.5, 'Renderer', 'Painters');

    subplot(3,1,1:2),hold on,
       plot(Time/60, NucTFEB, 'LineWidth', 2.50, 'Color', [.30,0.75,0.93])
       ylabel('Nuclear TFEB (%)')
       title('TFEB nuclear translocation')
       set(gca, 'XLim', [0,length(in)*15-15], 'YLim', [0,1], 'XTick',...
           0:100:length(in)*15, 'XTickLabel', {},'Box','off', 'YTick',...
           0:.2:1, 'YTickLabel', 0:20:100);

    subplot(3,1,3)
    plot([0:Ts:dimEXP-1], in, 'LineWidth', 2.5, 'Color', [.99 .55 .38]);
        ylabel('Input')
        xlabel('Time (min)'); 
        set(gca, 'XLim', [0,length(in)*15-15], 'XTick',  0:100:length(in)*15, 'XTickLabel', ...
            {0:100:length(in)*15}, 'YLim', [-.1 1.1], 'Box', 'off', 'YTick', ...
            [0 1], 'YTickLabel', {'RPMI','HBSS'});
        
    print(F,strcat(ExpName,'_',Mod),'-dpng')
    savefig(F,strcat(ExpName,'_',Mod,'.fig'))
    
    close
    
end