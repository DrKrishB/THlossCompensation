% Index of overlapping neurons that may be the same
clear P % P is the cell inside which merged neurons are stored
x = zeros(n,1) ; y = zeros(n,1) ; 
for j = 1:n
    ctd = S{j}.xyz ;
    x(j) = ctd(1) ;
    y(j) = ctd(2) ;
end
merge = zeros(n,1) ;
for j = 1:n
    if merge(j) == 0
        dst = zeros(n,1) ;
        for i = 1:n
            dst(i) = round(sqrt((x(j)-x(i))^2 + (y(j) - y(i))^2)) ;
        end
        idx = find(dst <= min_dist) ;
        nd = length(idx) ;
        if nd > 1
            for i = 1:nd
                if j ~= idx(i)
                    if round(corr(S{j}.Fraw,S{idx(i)}.Fraw),2) >= min_corr  
                        merge(idx(i)) = j ;                                
                    end
                end
            end
        end
    end
end
% merge neurons
k = 0; wk4 = []; wk5 = [];
for j = 1:n
    if ~isnan (merge(j))
        idx = find (merge == j) ; 
        if length(idx) == 0   
            k = k+1 ;
            P{k} = S{j} ;
        else % create a composite neuron from all overlaps 
            idx = sort(unique([j; idx])) ;    
            scs = []; id = []; cflt=[]; craw = [];
            for jk = 1:length(idx)
                id{jk} = S{idx(jk)}.ID ;
            end
            imb = zeros(d1,d2,length(idx)) ;
            for jk = 1:length(idx)
                % construct a mask            
                imb(:,:,jk) = S{idx(jk)}.img ;                 
                scs = [scs, S{idx(jk)}.score] ; % classifier scores
                cflt = cat(2,cflt,S{idx(jk)}.Cflt) ;
                craw = cat(2,craw,S{idx(jk)}.Craw) ;
                if jk > 1
                    merge(idx(jk)) = nan ;
                end
            end 
            a = imbilatfilt(max(imb,[],3)) ;                   
            bw = imfill(imbinarize(a),'holes') ;
            %bw = bwareafilt(bw,[minArea,20*minArea]) ;   
            s = regionprops(bw, a, 'All') ;
            if size(s,1) > 1      % choose the brightest & most circular  
                [I1,idx] = sort(cat(1,s.Circularity).*cat(1,s.MeanIntensity), 'descend') ;   
                bw = zeros(d1*d2,1) ;
                bw(s(idx(1)).PixelIdxList) = 1 ;
                bw = reshape(bw, [d1,d2]) ; 
                s = regionprops(bw, a, 'All') ;
            end
            if ~isempty(s)  
                k = k+1 ;
                P{k}.ID = cell2mat(id) ;
                P{k}.stdproj = mean(movstd(s.PixelIdxList)) ;
                P{k}.meanproj = mean(movavg(s.PixelIdxList)) ;
                P{k}.maxproj = mean(movmax(s.PixelIdxList)) ;
                P{k}.center = round(s.WeightedCentroid) ;
                P{k}.score = scs ;
                P{k}.Cflt = mean(cflt,2) ;
                P{k}.Craw = mean(craw,2) ;
                temp = zeros(m,1) ;
                for i = 1:m        
                    temp(i) = mean(mov(s.PixelIdxList,i)) ;
                end
                P{k}.Fraw = temp ;
                %imb = zeros(d1*d2,1) ;
                %a = reshape(a,[d1*d2,1]) ;
                %imb(s.PixelIdxList) = a(s.PixelIdxList) ;
                P{k}.img = sparse(a) ;
                P{k}.xyz = [pixel*s.WeightedCentroid,(np-1)*stepZ] ;
            else
                wk5 = [wk5 id] ;
            end
        end
    else
        wk4 = [wk4 S{j}.ID] ;
    end
end
writematrix(wk4', cat(2, outnam, ' rejected due to merging.csv'), 'WriteMode','append') ; 
if ~isempty(wk5)
    writematrix(cell2mat(wk5)', cat(2, outnam, ' rejected after merging.csv'), 'WriteMode','append') ;
end