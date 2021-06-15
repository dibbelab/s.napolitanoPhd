%%                           TRACKING CELLS.
%                            Sara Napolitano
%                               21/6/2019

function CELLS=trackingCELLS(CELLS,path,fields,i)

    MAX_NUCLEI_DIMENSION = 100;
        
    for j = 1:fields
        Cells = [CELLS{1,j},CELLS{2,j}];
        RemovedCells = [];
        idN = CELLS{3,j};
        
        cropfile = strcat('CropRect',num2str(j),'.mat');
        load(cropfile);
        
        fileGREEN = strcat(path,'t',num2str(i,'%.2d'),'xy',num2str(j,'%.2d'),'c2.tif');
        fileRED = strcat(path,'t',num2str(i,'%.2d'),'xy',num2str(j,'%.2d'),'c3.tif');
        fileMASK = strcat(path,'segmentation\t',num2str(i,'%.2d'),'xy',num2str(j,'%.2d'),'c1_m00_mask.png');
        
        s1=0; s2=0; s3=0;
        
        while (s2==0)||(s3==0)||(s1==0)
            s1 = ExistImage(fileMASK);
            s2 = ExistImage(fileGREEN);
            s3 = ExistImage(fileRED);
        end
        
        if s1==1&&s2==1&&s3==1
    
            Im_Green = imcrop(imread(fileGREEN),CropRect);
            Im_Red = imcrop(imread(fileRED),CropRect);
            Im_Mask = imcrop(imread(fileMASK),CropRect);
            
            % BACKGROUND REMOVAL
            background = imopen(Im_Green,strel('disk',40));
            Im_Green = imsubtract(Im_Green,background);
            
            background = imopen(Im_Red,strel('disk',40));
            Im_Red = imsubtract(Im_Red,background);
            
            % THE MASK IMAGE HAS TO BE LOGIC
            Mask = false(size(Im_Mask));
            Mask(Im_Mask ~= 0) = true;
            
            % NUCLEI SEGMENTATION
            stats = FindNuclei(Im_Red);
            statscell = regionprops(Mask,'Centroid','MajorAxisLength','MinorAxisLength','PixelIdxList');
            
            if ~isempty(stats)
            
                if i==1
                    
                    for k=1:length(stats)
                        POS = 1;
                        Cells(k).label = idN;
                        idN = idN+1;
                        
                        Cells(k).frame(POS) = i;
                        
                        BW = false(size(Im_Green));
                        BW(stats(k).PixelIdxList) = 1;
                        NucleusGreenFluo = regionprops(BW,Im_Green,'PixelValues');  % 'MeanIntensity',
                        
                        [CytoFluo,CellFluo,PixelIdxList,Area] = CytoFluoEst(stats(k),statscell,Im_Green,Cells(k));    %MeanCytoFluo,MeanCellFluo,
                        
                        clear BW
                        
                        Cells(k).Centroid(POS,:) = stats(k).Centroid; 
                        Cells(k).MajorAxisLength(POS) = stats(k).MajorAxisLength; 
                        Cells(k).MinorAxisLength(POS) = stats(k).MinorAxisLength; 
                        Cells(k).Orientation(POS) = stats(k).Orientation;
                        Cells(k).AreaNucleus(POS) = stats(k).Area;
                        Cells(k).CellArea(POS) = Area;
                        Cells(k).frame(POS) = i;
                        Cells(k).PixelIdxList{POS} = stats(k).PixelIdxList;
                        Cells(k).CellPixelIdxList{POS} = PixelIdxList;
                        Cells(k).GreenFluo(POS) = sum(NucleusGreenFluo.PixelValues);
                        Cells(k).CytoFluo(POS) = CytoFluo;
                        Cells(k).CellFluo(POS) = CellFluo;
                        Cells(k).RedFluo(POS) = sum(stats(k).PixelValues);
                        Cells(k).nucleusLABEL = 0;
                        Cells(k).originalCELL = Cells(k).label;
                        Cells(k).field = j;
                       
                        for q=1:length(Cells)
                            if norm([Cells(q).Centroid(end,1) Cells(q).Centroid(end,2)] - [Cells(k).Centroid(end,1) Cells(k).Centroid(end,2)]) < MAX_NUCLEI_DIMENSION
                                Cells(k).nucleusLABEL = 1;
                                Cells(k).originalCELL = Cells(q).label;
                            end
                        end
                        
                    end

                else 
                    
                    OldObj = extractOBJS(Cells);
                    NewObj = extractOBJS(stats);
                    [m,status1] = solveLP(OldObj,NewObj);

                    if status1==1
                        [Cells,RemovedCells,idN] = updateOBJS(Im_Green,m,Cells,RemovedCells,stats,statscell,idN,i,j);
                    end
                    
                end
                
            else
                
                ImSize = size(Im_Red);
                maskNUCLEI = false(ImSize);
                RemovedCells = [RemovedCells Cells];
                
            end

        end
        
        CELLS{1,j} = Cells;
        CELLS{2,j} = RemovedCells;
        CELLS{3,j} = idN;
        
    end
    
end