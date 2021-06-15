%%                             FIND NUCLEI
%                            Sara Napolitano
%                              22/07/2019

function [CytoFluo,CellFluo,PixelIdxList,CellArea] = CytoFluoEst(stats,statscell,Im_Green,Cell)

%     MeanCytoFluo = NaN;
%     MeanCellFluo = NaN;
    CytoFluo = NaN;
    CellFluo = NaN;
    CellArea = NaN;
    PixelIdxList = [];
    
    BWfaster = false(size(Im_Green));
    
    BWnuclei = false(size(Im_Green));
    BWnuclei(stats.PixelIdxList) = true;
    
    for index = 1:length(statscell)
        BWfaster = false(size(Im_Green));
    
        if norm(statscell(index).Centroid - stats.Centroid) <= statscell(index).MinorAxisLength/2
            BWfaster(statscell(index).PixelIdxList) = true;
            BWcell = BWfaster | (BWnuclei);
            BWcyto = BWcell & ~(BWnuclei);
            
            clear BWfaster s
            
            statsCell = regionprops(BWcell,Im_Green,'Area','PixelIdxList','PixelValues'); % ,'MeanIntensity'
            statsCyto = regionprops(BWcyto,Im_Green,'Area','PixelIdxList','PixelValues');   % ,'MeanIntensity'
            
            PixelValues = statsCell.PixelValues;
            CellFluo = sum(PixelValues);
            CellArea = statsCell.Area;
            PixelIdxList = statsCell.PixelIdxList;
            
           if length(statsCyto) == 1
                CytoFluo = sum(statsCyto.PixelValues);
            else
                Area = 0;
                FluoNotNorm = 0;
                for index_cyto = 1:length(statsCyto)
                    Area = Area + statsCyto(index_cyto).Area;
                    FluoNotNorm = FluoNotNorm + sum(statsCyto(index_cyto).PixelValues);
                end
                CytoFluo = FluoNotNorm;
            end
        end
    end
    
    if isnan(CellFluo) && length(Cell.frame) > 1  
       
        BWfaster(Cell.CellPixelIdxList{end}) = true;
        
        BWcell = BWfaster | (BWnuclei);
        BWcyto = BWcell & ~(BWnuclei);
            
        statsCell = regionprops(BWcell,Im_Green,'Area','PixelIdxList','PixelValues'); % ,'MeanIntensity'
        statsCyto = regionprops(BWcyto,Im_Green,'Area','PixelIdxList','PixelValues');   % ,'MeanIntensity'

        if ~isempty(statsCell)
    
            PixelValues = statsCell.PixelValues;
            CellFluo = sum(PixelValues);
            CellArea = statsCell.Area;
            PixelIdxList = statsCell.PixelIdxList;

            if length(statsCyto) == 1
                CytoFluo = sum(statsCyto.PixelValues);
            else
                Area = 0;
                FluoNotNorm = 0;
                for index_cyto = 1:length(statsCyto)
                    Area = Area + statsCyto(index_cyto).Area;
                    FluoNotNorm = FluoNotNorm + sum(statsCyto(index_cyto).PixelValues);
                end
                CytoFluo = FluoNotNorm;
            end
        end
    end
    
end