%% 

function [] = makeFig(expNAME,inputLABEL,inputTICK)

    if isfile(strcat(pwd,'\Workspaces\',expNAME,'_toClean.mat'))
        
        load(strcat(pwd,'\Workspaces\',expNAME,'_toClean.mat'));

        % -----------------------------------------------------------------
        % to clear the data-set
        cells(toREMOVE) = [];
        DeadTimes(toREMOVE) = [];
        DivisionTime(:,toREMOVE) = [];
        
        if length(in) >= max(DeadTimes)
            in = in(1:max(DeadTimes));
        else
            in = [in;in(end).*ones(max(DeadTimes)-length(in),1)];
        end

        clear toREMOVE 

        fluoMA = NaN(length(cells),max(DeadTimes));
        cytofluoMA = NaN(length(cells),max(DeadTimes));

        for i = 1:length(DeadTimes)

            ind = find(cells(i).frame == DeadTimes(i));

            if isempty(ind)
                ind = length(cells(i).frame);
            end

            if ~isnan(DivisionTime(1,i))
                ind1 = find(cells(i).frame == DivisionTime(1,i));
                ind2 = find(cells(i).frame == DivisionTime(2,i));

                if isempty(ind1)
                    ind1 = 1;
                end

                temp = cells(i).GreenFluo(1:ind);
                temp(ind1:ind2) = NaN;

                TEMP = cells(i).CellFluo(1:ind);
                TEMP(ind1:ind2) = NaN;

                Temp = cells(i).CytoFluo(1:ind);
                Temp(ind1:ind2) = NaN;

                if ~isnan(DivisionTime(3,i))
                    ind3 = find(cells(i).frame == DivisionTime(3,i));
                    ind4 = find(cells(i).frame == DivisionTime(4,i));

                    temp(ind3:ind4) = NaN;
                    TEMP(ind3:ind4) = NaN;
                    Temp(ind3:ind4) = NaN;
                end
                
                if size(DivisionTime,1) == 6
                    if ~isnan(DivisionTime(5,i))
                        ind5 = find(cells(i).frame == DivisionTime(5,i));
                        ind6 = find(cells(i).frame == DivisionTime(6,i));

                        temp(ind5:ind6) = NaN;
                        TEMP(ind5:ind6) = NaN;
                        Temp(ind5:ind6) = NaN;
                    end
                end

                M = movmean(TEMP,20,'omitnan');

                toSaveMA = temp./M;
                cytoMA = Temp./M;

            else

                M = movmean(cells(i).CellFluo(1:ind),20,'omitnan');

                toSaveMA = cells(i).GreenFluo(1:ind)./M;
                cytoMA = cells(i).CytoFluo(1:ind)./M;

            end

            if max(toSaveMA) > 1
                toSaveMA = toSaveMA-(max(toSaveMA)-1);
                toSaveMA(toSaveMA<0) = 0;
            end

            fluoMA(i,cells(i).frame(1:ind)) = toSaveMA;
            cytofluoMA(i,cells(i).frame(1:ind)) = cytoMA;

        end

        clear ind1 ind2 ind3 ind4 ind5 ind6 toSaveMA M ind DeadTimes cells 
        clear temp Temp TEMP cytoMA i DivisionTime

        save(strcat(pwd,'\Workspaces\',expNAME,'_Final.mat'));


        % -----------------------------------------------------------------

        F = figure('Position', [1 1 720 720], 'DefaultAxesFontSize', 12, ...
                'DefaultAxesLineWidth', 2.5, 'Renderer', 'Painters');

        subplot(5,1,1:2),hold on,
        for k = 1:size(fluoMA,1)
            plot([0:length(in)-1]*15,fluoMA(k,:),'Color',[0.51 0.80 0.35],...
                'LineWidth', .5);
        end
        hold on, plot([0:length(in)-1]*15,mean(fluoMA,'omitnan'),...
            'LineWidth', 2.5, 'Color', [0.18,0.46,0.15]);
        ylabel('Nuclear TFEB (%)')
        title('TFEB nuclear translocation')
        set(gca, 'XLim', [0,length(in)*15-15], 'YLim', [0,1], 'XTick',...
            0:100:length(in)*15, 'XTickLabel', {},'Box','off', 'YTick',...
            0:.2:1, 'YTickLabel', 0:20:100);

        subplot(5,1,3:4),hold on,
        for k = 1:size(cytofluoMA,1)
            plot([0:length(in)-1]*15,cytofluoMA(k,:),'Color',[0.51 0.80 0.35],...
                'LineWidth', .5);
        end
        hold on, plot([0:length(in)-1]*15,mean(cytofluoMA,'omitnan'),...
            'LineWidth', 2.5, 'Color', [0.18,0.46,0.15]);
        ylabel('Cytosolic TFEB (%)')
        set(gca, 'XLim', [0,length(in)*15-15], 'YLim', [0,1], 'XTick',...
            0:100:length(in)*15, 'XTickLabel', {},'Box','off', 'YTick',...
            0:.2:1, 'YTickLabel', 0:20:100);

        subplot(5,1,5), plot([0:length(in)-1]*15,in,...
            'LineWidth', 2.5, 'Color', [.99 .55 .38]);
        ylabel(inputLABEL)
        xlabel('Time (min)'); 
        set(gca, 'XLim', [0,length(in)*15-15], 'XTick',  0:100:length(in)*15, 'XTickLabel', ...
            {0:100:length(in)*15}, 'YLim', [-.1 1.1], 'Box', 'off', 'YTick', ...
            [0 1], 'YTickLabel', inputTICK);

        print(F,strcat(pwd,'/Figures/',expNAME),'-dpng')
        savefig(F,strcat(pwd,'/Figures/',expNAME,'.fig'))

        close

        clear F in cytofluoMA fluoMA

    else
        
        load(strcat(pwd,'\Workspaces\',expNAME,'_Final.mat'));
        
        F = figure('Position', [1 1 720 435], 'DefaultAxesFontSize', 12, ...
                'DefaultAxesLineWidth', 2.5, 'Renderer', 'Painters');

        subplot(3,1,1:2),hold on,
        for k = 1:size(fluoMA,1)
            plot([0:length(in)-1]*15,fluoMA(k,:),'Color',[0.51 0.80 0.35],...
                'LineWidth', .5);
        end
        hold on, plot([0:length(in)-1]*15,mean(fluoMA,'omitnan'),...
            'LineWidth', 2.5, 'Color', [0.18,0.46,0.15]);
        ylabel('Nuclear TFEB (%)')
        title('TFEB nuclear translocation')
        set(gca, 'XLim', [0,length(in)*15-15], 'YLim', [0,1], 'XTick',...
            0:100:length(in)*15, 'XTickLabel', {},'Box','off', 'YTick',...
            0:.2:1, 'YTickLabel', 0:20:100);

        subplot(3,1,3), plot([0:length(in)-1]*15,in,...
            'LineWidth', 2.5, 'Color', [.99 .55 .38]);
        ylabel(inputLABEL);
        xlabel('Time (min)'); 
        set(gca, 'XLim', [0,length(in)*15-15], 'XTick',  0:100:length(in)*15, 'XTickLabel', ...
            {0:100:length(in)*15}, 'YLim', [-.1 1.1], 'Box', 'off', 'YTick', ...
            [0 1], 'YTickLabel', inputTICK);

        print(F,strcat(pwd,'/Figures/',expNAME),'-dpng')
        savefig(F,strcat(pwd,'/Figures/',expNAME,'.fig'))

        close

        clear F in cytofluoMA fluoMA
    
    end
    
end