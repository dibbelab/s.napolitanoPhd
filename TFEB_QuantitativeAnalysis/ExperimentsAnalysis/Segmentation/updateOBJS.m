%%                            UPDATE OBJETS
%                            Sara Napolitano
%                               21/6/2019

% -------------------------------------------------------------------------
% è una funzione che aggiorna le cellule tra 2 frame
% -------------------------------------------------------------------------

function [Cells,removedCELLS,idN] = updateOBJS(Im_Green,m,Cells,removedCELLS,stats,statscell,idN,frame,field)
    
    MAX_NUCLEI_DIMENSION=100;
    M = size(m,1);
    updatedINDEXES = [];
    
    for k = 1:M
        INDEX = find(m(k,:) > 0);
    
        if length(INDEX) == 1
            
            Z = INDEX(1);
            POS = size(Cells(k).Centroid,1)+1;
           
            Cells(k).frame(POS) = frame;

            BW = false(size(Im_Green));
            BW(stats(Z).PixelIdxList) = 1;
            NucleusGreenFluo = regionprops(BW,Im_Green,'MeanIntensity','PixelValues');
            
            [CytoFluo,CellFluo,PixelIdxList,Area] = CytoFluoEst(stats(Z),statscell,Im_Green,Cells(k));
            
            clear BW
            
            Cells(k).Centroid(POS,:) = stats(Z).Centroid; 
            Cells(k).MajorAxisLength(POS) = stats(Z).MajorAxisLength; 
            Cells(k).MinorAxisLength(POS) = stats(Z).MinorAxisLength; 
            Cells(k).Orientation(POS) = stats(Z).Orientation;
            Cells(k).AreaNucleus(POS) = stats(Z).Area;
            Cells(k).CellArea(POS) = Area;
            Cells(k).PixelIdxList{POS} = stats(Z).PixelIdxList;
            Cells(k).CellPixelIdxList{POS} = PixelIdxList;
            Cells(k).GreenFluo(POS) = sum(NucleusGreenFluo.PixelValues);
            Cells(k).CytoFluo(POS) = CytoFluo;
            Cells(k).CellFluo(POS) = CellFluo;
            Cells(k).RedFluo(POS) = sum(stats(Z).PixelValues);
            
            updatedINDEXES(end+1) = k;
        
        elseif length(INDEX) > 1
            
            NORMS = [];
        
            for h = 1:length(INDEX)
                NORMS(h) = norm([Cells(k).Centroid(end,1) Cells(k).Centroid(end,2)] - stats(INDEX(h)).Centroid(:));
            end

            [~,IDS] = sort(NORMS);
            sortedINDEXES = INDEX(IDS);
        
            POS = size(Cells(k).Centroid,1)+1;
            
            Cells(k).frame(POS) = frame;
            
            BW = false(size(Im_Green));
            BW(stats(sortedINDEXES(1)).PixelIdxList) = 1;
            NucleusGreenFluo = regionprops(BW,Im_Green,'MeanIntensity','PixelValues');
            
            [CytoFluo,CellFluo,PixelIdxList,Area] = CytoFluoEst(stats(sortedINDEXES(1)),statscell,Im_Green,Cells(k));
            
            clear BW

            Cells(k).Centroid(POS,:) = stats(sortedINDEXES(1)).Centroid; 
            Cells(k).MajorAxisLength(POS) = stats(sortedINDEXES(1)).MajorAxisLength; 
            Cells(k).MinorAxisLength(POS) = stats(sortedINDEXES(1)).MinorAxisLength; 
            Cells(k).Orientation(POS) = stats(sortedINDEXES(1)).Orientation;
            Cells(k).AreaNucleus(POS) = stats(sortedINDEXES(1)).Area;
            Cells(k).CellArea(POS) = Area;
            Cells(k).PixelIdxList{POS} = stats(sortedINDEXES(1)).PixelIdxList;
            Cells(k).CellPixelIdxList{POS} = PixelIdxList;
            Cells(k).GreenFluo(POS) = sum(NucleusGreenFluo.PixelValues);
            Cells(k).CytoFluo(POS) = CytoFluo;
            Cells(k).CellFluo(POS) = CellFluo;
            Cells(k).RedFluo(POS) = sum(stats(sortedINDEXES(1)).PixelValues);
                        
            updatedINDEXES(end+1) = k;
            
            for D = 2:length(sortedINDEXES)

                S = length(Cells);
                POS = 1;
                
                if S+1>length(Cells)
                    Cells(S+1).label = idN;
                    idN = idN+1;
                end
                
                Cells(S+1).field = field;
                
                Cells(S+1).frame(POS) = frame;
                
                BW = false(size(Im_Green));
                BW(stats(sortedINDEXES(D)).PixelIdxList) = 1;
                NucleusGreenFluo = regionprops(BW,Im_Green,'MeanIntensity','PixelValues');
                
                [CytoFluo,CellFluo,PixelIdxList,Area] = CytoFluoEst(stats(sortedINDEXES(D)),statscell,Im_Green,Cells(S+1));

                clear BW
                
                Cells(S+1).Centroid(POS,:) = stats(sortedINDEXES(D)).Centroid; 
                Cells(S+1).MajorAxisLength(POS) = stats(sortedINDEXES(D)).MajorAxisLength; 
                Cells(S+1).MinorAxisLength(POS) = stats(sortedINDEXES(D)).MinorAxisLength; 
                Cells(S+1).Orientation(POS) = stats(sortedINDEXES(D)).Orientation;
                Cells(S+1).AreaNucleus(POS) = stats(sortedINDEXES(D)).Area;
                Cells(S+1).CellArea(POS) = Area;
                Cells(S+1).PixelIdxList{POS} = stats(sortedINDEXES(D)).PixelIdxList;
                Cells(S+1).CellPixelIdxList{POS} = PixelIdxList;
                Cells(S+1).GreenFluo(POS) = sum(NucleusGreenFluo.PixelValues);
                Cells(S+1).CytoFluo(POS) = CytoFluo;
                Cells(S+1).CellFluo(POS) = CellFluo;
                Cells(S+1).RedFluo(POS) = sum(stats(sortedINDEXES(D)).PixelValues);
                
                if norm([Cells(k).Centroid(end,1) Cells(k).Centroid(end,2)] - [Cells(S+1).Centroid(end,1) Cells(S+1).Centroid(end,2)]) < MAX_NUCLEI_DIMENSION
                    Cells(S+1).nucleusLABEL=1;
                    Cells(S+1).originalCELL=Cells(k).label;
                end

                updatedINDEXES(end+1) = S+1;
            
            end
            
        end
        
    end
        
    toREMOVE = [];

    for h = 1:length(Cells)
    
        if isempty(find(updatedINDEXES==h,1))
        
            toREMOVE = [toREMOVE h];
    
        end
    
    end


    for h = 1:length(toREMOVE)
    
        removedCELLS = [removedCELLS Cells(toREMOVE(h))];
    
    end


    Cells(toREMOVE) = [];
    
end