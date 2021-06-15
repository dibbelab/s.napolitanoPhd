%% simulator.m
%%% JUNE, 2021

function [Time,NucTFEB] = simulator(Mod,InputLev,InitCond,dimExp)
    
    global u
    Time = [];
    xsave = [];
    
    switch Mod
        
        case 'OpenLoop'
            sel_ODEs = @OpenLoop_Fun;    %@(t,x) OpenLoop_Fun(t,x,u);
        case 'Feedback'
            sel_ODEs = @Feedback_Fun;    %@(t,x) Feedback_Fun(t,x,u);
            
    end
    
    % Set tSPAN
    tSPAN = (0:180:dimExp)*60;  % in s
    x0 = InitCond;
    
    % integration
    for ITR = 1:length(tSPAN)-1
        ts_i = tSPAN(ITR);
        ts_f = tSPAN(ITR+1);
        
        if ITR == 7
            u = InputLev(2);
        elseif mod(ITR,2)
            u = InputLev(1);
        else
            u = InputLev(2);
        end
        
        [ts,xs] = ode45(sel_ODEs,[ts_i,ts_f],x0);
        
        xsave = [xsave;xs(2:end,:)];
        Time = [Time;ts(2:end)];
        x0 = xs(end,:);
    end
    
    NucTFEB = xsave(:,2)+ xsave(:,3);

end