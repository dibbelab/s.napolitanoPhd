%%                             FIND NUCLEI
%                            Sara Napolitano
%                              22/07/2019

function stats = FindNuclei(Im_Red)
    
    MIN_DIMENSION = 10;
    AreaMean = 1400;
    s1 = []; s2 =[]; s = [];
    toREMOVE = [];
    
    Im = imadjust(Im_Red);
    BW = imbinarize(Im);
    BW = imfill(BW,'holes');

    s1 = regionprops(BW,Im_Red,{'Centroid','MajorAxisLength','MinorAxisLength',...
        'Orientation','Area','PixelIdxList','MeanIntensity','PixelValues'});
    
    for z = 1:length(s1)
        K = round(s1(z).Area/AreaMean);
        
        if K > 1
            toREMOVE = [toREMOVE; z];
            
            bw = false(size(BW));
            bw(s1(z).PixelIdxList) = true;
            
            % FROM PAPER "Separating Touching Cells using Pixel Replicated Elliptical Shape Models"
            % For more than a million replicate points use sampling at 10% instead.
            numPoints = sum(round(bwdist(~bw)));
            if ( numPoints < 1e6  )
                ptsReplicated = PixelReplicate(bw);
            else
                ptsReplicated = PixelReplicateSample(bw,0.1);
            end
 
            objPR = fitgmdist(ptsReplicated, K, 'replicates',5);
            if ~objPR.Converged
                objPR = fitgmdist(ptsReplicated, K, 'replicates',5);
            end
            
            % clustering
%             [x,y] = find(bw);
%             idx = objPR.cluster([y,x]);
            for kk = 1:K
                [~,mask_temp] = genEllipse(objPR.mu(kk,:),objPR.Sigma(:,:,kk),bw,sqrt(20/3));
                stemp = regionprops(mask_temp,Im_Red,{'Centroid','MajorAxisLength',...
                    'MinorAxisLength','Orientation','Area','PixelIdxList',...
                    'MeanIntensity','PixelValues'});
                s2 = [s2;stemp];
                clear mask_temp stemp
            end
        end
    end
    
    s1(toREMOVE) = [];
    
    s = [s1;s2];

    stats = [];

    for q = 1:length(s)
        
        toADD = true;
        
        if s(q).MinorAxisLength < MIN_DIMENSION
            toADD = false;
        end

        for h = 1:size(stats,1)
            
            if norm(stats(h).Centroid-s(q).Centroid) < MIN_DIMENSION
                toADD = false;
            end
            
        end

        if toADD
            stats = [stats;s(q)];
        end
        
    end
    
end