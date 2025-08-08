%% Remove false positive neurons & events based on activity
% normalize between 0 & 1 
B = zeros(n,m) ; C = B ;
for j = 1:n
    B(j,:) = Q{j}.Fraw ;
    C(j,:) = Q{j}.Cflt ;
end
B = mat2gray(B) ;
C = mat2gray(C) ;
baseline = zeros(n,1) ;
D = zeros(n,m) ;
% estimate baseline for each neuron
for j = 1:n
    tmp = B(j,:) ;
    tmpF = smoothdata(tmp, 'movmean', round(120*fps)) ;
    tmp = tmp - tmpF ;
    t = tmp ;
    t = t(t>0) ;        
    baseline(j) = rms(t(t < median(t))) ;
end
minh = rms(baseline) ;  
k = 0; 
for j = 1:n
    tmp = B(j,:)' ;
    tmpF = smoothdata(tmp, 'movmean', round(120*fps)) ;
    tmp = tmp - (tmpF+minh) ;
    Q{j}.dfF = tmp ;
    tms = sort(tmp,'descend') ;
    %if mean(tms(1:20)) > 2* max(minh,rms(tmp(tmp>0)))
        temp = zeros(m,1) ;
        for i = 1:m
            if tmp(i) > 0                       
               temp(i) = tmp(i) ;
            end
        end
        % remove single-frame events
        for i = 2:m-1   
            if temp(i) > 0 && temp(i-1) == 0 && temp(i+1) == 0 
                temp(i) = 0;
            end
        end
        % remove dual-frame events
        for i = 3:m-2   
            if temp(i) > 0 && temp(i+1) > 0
                if temp(i-1) == 0 && temp(i+2) == 0 
                    temp(i) = 0;
                    temp(i+1) = 0;
                end
            elseif temp(i) > 0 && temp(i-1) > 0
                if temp(i-2) == 0 && temp(i+1) == 0 
                    temp(i) = 0;
                    temp(i-1) = 0;
                end
            end
        end
        Q{j}.Fflt = temp ;
        if length(find(temp>0)) > 1
            k = k+1 ;
            R{k} = Q{j} ;                
        end
    %end        
end