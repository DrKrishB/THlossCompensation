k = 0 ; 
bd = round(0.02*mean(d1,d2)) ; % distance from frame boundaries inside which the neuron is to be rejected
mov = reshape(mov, [d1*d2,m]) ;
%maxproj_thresh = rms(movmax(movmax<mean(movmax))) ;
%otsu_thresh = (0.5+0.015*np)*graythresh(ref) ; 
maxproj_thresh = thf*YenThreshKB(ref,4096)/4096 ;
maxproj_thresh = min(0.2,maxproj_thresh) ;
movstd = reshape(movstd, [d1*d2,1]) ;  % 1D Standard deviation projection 
movavg = reshape(movavg, [d1*d2,1]) ;  % 1D Average projection
movmax = reshape(movmax, [d1*d2,1]) ;  % 1D Maximal projection
movmin = reshape(movmin, [d1*d2,1]) ;  % 1D Minimum projection
ref = reshape(ref, [d1*d2,1]) ;  % 1D Max-Minimum projection

wk1 = [] ; %weak neurons
wk2 = [] ; %neurons rejected due to flatness
wk3 = [] ; %neurons rejected due to morphology
for j = 1:length(stat)   
    X = stat{j}.xpix + 1 ; % 1 is added as in python index starts with 0
    Y = stat{j}.ypix + 1 ;
    I = stat{j}.lam ;   % pixel values
    p = sub2ind([d1,d2],Y,X) ; 
    % reject neuron if mean intensity is less than maxproj_thresh
    if mean(ref(p)) >= maxproj_thresh
        %bw = zeros(d1*d2,1) ;
        %bw(p) = 1 ;
        %bw = imfill(imdilate(reshape(bw,[d1,d2]),strel('disk',1)),'holes') ;
        %p = find(bw>0) ;
        a = zeros(d1*d2,1) ;
        for i = 1:length(p)
            %a(p(i)) = I(i) ; %Make neuron using suite2p's spatial component
            %a(p(i)) = ref(p(i)) ; %Make neuron using max-min projection 
            a(p(i)) = mean([I(i); 5*ref(p(i))]) ;
        end
        a = imbilatfilt(imdilate(reshape(a, [d1,d2]),strel('disk',1))) ; 
        a1 = mat2gray(a) ;
        bw = imbinarize(a1,0.5) ;   
        bw1 = imbinarize(a1,0.4) ;
        bw2 = imbinarize(a1,0.6) ;
        % For a neuron that is real and robust, the size of non-zero pixels in
        % bw1 & bw2 should not vary by more than 50% of that in bw
        pxc = [length(find(bw1>0)),length(find(bw>0)),length(find(bw2>0))] ;
        if mean(abs(diff(pxc))) <= 1+ceil(0.5*length(find(bw>0)))
            % properties of the neuron
            bw = bwareafilt(bw,[minArea,15*minArea]) ;            
            s = regionprops(bw, a, 'All') ;
            if size(s,1) > 1      
                % choose the brightest & most circular  
                %[I1,idx] = sort(cat(1,s.Circularity).*cat(1,s.MeanIntensity), 'descend') ; 
                % choose the brightest & largest  
                [I1,idx] = sort(cat(1,s.Area).*cat(1,s.MeanIntensity), 'descend') ;
                bw = zeros(d1*d2,1) ;
                bw(s(idx(1)).PixelIdxList) = 1 ;
                bw = reshape(bw, [d1,d2]) ; 
                s = regionprops(bw, a, 'All') ;
            end
             
            if ~isempty(s)  
                bw = zeros(d1*d2,1) ;
                cn = round(s.WeightedCentroid) ;                 
                %if ismember(sub2ind([d1,d2],cn(2),cn(1)), mps)
                if cn(1) < d2-bd && cn(2) < d1-bd && cn(1) > bd && cn(2) > bd
                    if s.MajorAxisLength > minLength && s.MinorAxisLength > minWidth 
                        if round(s.Eccentricity,2) <= maxEcn                    
                            bw(s.PixelIdxList) = 1 ;                     
                        end
                    end
                end
                %end 
                bw = reshape(bw,[d1,d2]) ;
                bw = imbinarize(bw) ;
                s = regionprops(bw, a, 'All') ;
            end            

            if ~isempty(s)        
                k = k+1 ;
                S{k}.ID = j ;
                S{k}.stdproj = mean(movstd(s.PixelIdxList)) ;
                S{k}.meanproj = mean(movavg(s.PixelIdxList)) ;
                S{k}.maxproj = mean(movmax(s.PixelIdxList)) ;
                S{k}.center = round(s.WeightedCentroid) ;
                S{k}.score = round(iscell(j,2),2) ;
                if length(spks(j,:)) == 1401
                    S{k}.Cflt = spks(j,2:1401)' ;
                    S{k}.Craw = F(j,2:1401)' ;      
                else
                    S{k}.Cflt = spks(j,:)' ;
                    S{k}.Craw = F(j,:)' ; 
                end
                temp = zeros(m,1) ;
                for i = 1:m        
                    temp(i) = mean(mov(s.PixelIdxList,i)) ;
                end
                S{k}.Fraw = temp ;
                imb = zeros(d1*d2,1) ;
                a = reshape(a,[d1*d2,1]) ;
                imb(s.PixelIdxList) = a(s.PixelIdxList) ;
                S{k}.img = sparse(reshape(imb,[d1,d2])) ;
                S{k}.xyz = round([pixel*s.WeightedCentroid,(np-1)*stepZ]) ;
            else
                wk3 = [wk3 j] ;
            end
        else
            wk2 = [wk2 j] ;
        end
    else
        wk1 = [wk1 j] ;
    end
end  
csvwrite(cat(2, outnam, ' rejected due to mean intensity.csv'), wk1') ;     
csvwrite(cat(2, outnam, ' rejected due to flatness.csv'), wk2') ;    
csvwrite(cat(2, outnam, ' rejected due to morphology.csv'), wk3') ;