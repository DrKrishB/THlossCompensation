%% Author: Dr. Krishnashish Bose (Guo Lab, UCSF)
%% Created on: 22-May-2020
%% Last modified on: 16-Jul-2024 
%% To run on Unix, gstreamer needs to be installed
%% Status: Completed 
function [] = track240fps(fn)
% clear; close all; clc
warning('off')
ext = '.MOV' ; % MOV supported in Windows 11 only, otherwise converted using ffmpeg
juvenile = 0  % 1 => juvenile larvae. 0 => <17 dpf larvae
pigmented = 0  % 1 => Pigmented larvae. 0 => Non-pigmented larvae
saveavi = 1  % 1 => creates 30s movie clips (takes more space and time)
savemat = 1  % 1 => saves a mat file for each arena (takes more time) 
border = 0     % 1 => objects touching the walls are rejected. Not recommended.
videoquality = 70 ; % quality of saved tracking video files
j1 = 1; j2 = 0 ;
px = 0.14 ; % pixel size in mm
arenasize = 34.5 ; %mm 
bdl = 4.5 ; % Body length (mm) of 7-8 dpf zebrafish larva
bgopt = 2 ; % background estimated from (1) first & last frBG frames (fast) 
            % or (2) every second

% user-specified cut-offs in mm
if juvenile == 1
    minWidth = 0.5 ;  % expected minimum width of a single eye in pixels
    maxWidth = 3 ;  % expected maximum width of the cluster containing both eyes
    minLength = 0.5 ; % expected maximum length of the cluster containing entire body
    if pigmented == 1    
        maxLength = 20 ; % expected maximum length of the cluster containing containing entire body
        minArea = 35 ;   % minimum pixels that might represent the fish    
        maxArea = 10*minArea ;
    else
        maxLength = 15 ;
        minArea = 15 ;   
        maxArea = 10*minArea ;
    end
else
    minWidth = 0.3 ;  % minimum width of a single eye in pixels
    maxWidth = 1 ;    % maximum width of the cluster containing both eyes
    minLength = 0.4 ; % maximum length of the cluster containing entire body
    if pigmented == 1    
        maxLength = 4 ; % maximum length of the cluster containing entire body
        minArea = 15 ;   % minimum pixels that might represent the fish
        maxArea = 150 ;   % maximum pixels that might represent the fish  
    else
        maxLength = 2 ;
        minArea = 9 ;
        maxArea = 90 ;
    end
end
% delete(gcp('nocreate'))
% parpool(6) ;
start = datetime   
file = dir(cat(2,'*',ext)) ;
% for fn = 1:length(file)    
tic
% if length(file) > 1
%     fprintf('%d MOV files found in %s \n', length(file), pwd)
%     fn = input('Which file?  ');  %file no.
% else
%     fn = 1 ;
% end
fname = file(fn).name ;    
% m = 15*60  % frames to be analyzed
% t1 = input('\nStart frame (mm.ss) : ','s') ;
% t2 = input('\nEnd frame (mm.ss) : ','s') ;  
% nam = cat(2,erase(fname,ext),' ',t1,'-',t2)
nam = erase(fname,ext) ;  

% ts = strsplit(t1,'.') ;
% ts1 = str2num(ts{1}) ;
% ts2 = str2num(ts{2}) ;
% m1 = 60*ts1 + ts2 ;
m1 = 0 ;
% ts = strsplit(t2,'.') ;
% ts1 = str2num(ts{1}) ;
% ts2 = str2num(ts{2}) ;
% m2 = 60*ts1 + ts2 ;
m2 = 0 ;
%matlab.video.read.UseHardwareAcceleration('off') ;
V = VideoReader(fname) ;
D1 = V.Width; D2 = V.Height;
fps = V.FrameRate ;
% fps = round(fps) 
frame_buffer = round(fps*30) ; %only used for saving video clips
frBG = round(fps*60) ; % 1 minute for background estimation (limited by RAM)
    
ml = V.NumFrames ;
m = round(900*fps) ; % I need to analyze 900 seconds

if m > ml || V.Duration < 900
    fprintf('\n%s does not have 15 minutes of recording', nam) ;
    % break ;
end
if m2 == 0
    if m1 == 0 
        m1 = floor((ml-m)/2) ;
    end
    m2 = m1 + m - 1 ;
end

fprintf('Frames %d to %d will be analyzed in %s recorded at %1.4f fps\n', m1,m2,nam,fps)
fps = round(fps) ;
%px = min(plate_size)*25.4/min(D1,D2) ; % pixel size in mm
% Background 
bgfile = cat(2,nam,' background.csv') ;
if isfile(bgfile)
    bgim = readmatrix(bgfile) ;
else
    fprintf('.........Estimating background .......\n')     
    if bgopt == 1 
        bgim = mean(read(V,m1),3) ;
        % first 1 minute
        while hasFrame(V) && i < frBG+fps            
            bgim = max(cat(3,bgim,mean(readFrame(V),3)),[],3) ;
            i = i+1 ;
        end 
        i = fps+m2-frBG ;
        img = read(V,i) ;
        % last 1 minute
        while hasFrame(V) && i < fps+m2
            bgim = max(cat(3,bgim,mean(readFrame(V),3)),[],3) ;            
            i = i+1 ;
        end 
    else % get the frame every 10 seconds
        bgim = mean(read(V,fps),3) ;
        for i = m1:10*fps:m2
            bgim = max(cat(3,bgim,mean(read(V,i),3)),[],3) ;
        end
    end
    bgim = bgim(150:1800,10:1070) ;
    writematrix(bgim,cat(2,nam,' background.csv'))
    %imwrite(uint8(bgim),cat(2,nam,' background.tif'))
end

% Automatic arena detection    
%     figure('Position',[1 1 2000 1500]) 
%     subplot(131)
%     imshow(bgim,[]); colorbar
 
 IMG = bgim ;
 IMG(IMG > 150) = 255 ;
%     subplot(132)
%     imshow(IMG, []); colorbar 
IMC = mat2gray(imbilatfilt(imcomplement(IMG))) ;
BW = imcomplement(imbinarize(IMC)) ; 
BW = bwpropfilt(BW, 'MinorAxisLength', round([30/px,50/px])) ;
BW = bwpropfilt(BW, 'MajorAxisLength', round([30/px,50/px])) ;
%     subplot(133)
%     imshow(bgim, []); colorbar
%     hold on ;
%     visboundaries(BW,'LineWidth',2,'EnhanceVisibility',false)
%     hold off ;
BB = regionprops(BW,'BoundingBox') ;
BB = cat(1,BB.BoundingBox) ;    
CC = regionprops(BW,'Centroid') ;
CC = cat(1,CC.Centroid) ;
n = size(CC,1) ;    % no. of wells in the arena
if j2 == 0
    j2 = n ;
end
CC = round(cat(2,CC,BB)) ;
[I,J1] = sortrows(round(0.01*CC),[2 1]) ;
q = I(:,2) ;
qid = find(diff(q)==1) ;
for ii = 1:length(qid)
    if mod(qid(ii),3) > 0
        q(qid(ii)) = q(qid(ii)) + 1 ;
    end
end
I(:,2) = q ;
%order = [1 2 7 8 3 4 9 10 5 6 11 12 13 14 19 20 15 16 21 22 17 18 23 24] ;
[I,J2] = sortrows(I,[2 1]) ;
CC = CC(J1(J2),:) ;  

%CC = CC(order,:) ;
imshow(bgim, []); colorbar
hold on ;
for j = j1:j2
    d1 = round(px*CC(j,5)) ; %mm size of square well
    d2 = round(px*CC(j,6)) ; %mm size of square well
    nd1 = d1-arenasize ;
    if nd1 > 0
        CC(j,3) = CC(j,3) + round(nd1/2) + 1 ;
        CC(j,5) = round(arenasize/px) - 1 ;
        CC(j,1) = CC(j,3) + CC(j,5)/2 - 1 ;
    end
    nd2 = d2-arenasize ;
    if nd2 > 0
        CC(j,4) = CC(j,4) + round(nd2/2) + 1 ;
        CC(j,6) = round(arenasize/px) - 1 ;
        CC(j,2) = CC(j,4) + CC(j,6)/2 - 1  ;
    end
    d1 = round(px*CC(j,5)) ; %mm size of square well
    d2 = round(px*CC(j,6)) ; %mm size of square well
    x0 = CC(j,3) ; y0 = CC(j,4) ;
    xd = CC(j,5) ; yd = CC(j,6) ;
    
    rectangle('Position', [x0 y0 xd yd], 'EdgeColor',[0.3,1,1]) ;
    x0 = x0+round(1.5*bdl/px) ; xd = xd-2*round(1.5*bdl/px)+1 ;
    y0 = y0+round(1.5*bdl/px) ; yd = yd-2*round(1.5*bdl/px)+1 ;
    rectangle('Position', [x0 y0 xd yd], 'EdgeColor',[0.4,0.3,0.8],'LineStyle','--') ;
    if j <10
        text(CC(j,1),CC(j,2),num2str(j),'Color','r') ;
    else
        text(CC(j,1)-10,CC(j,2),num2str(j),'Color','r') ;
    end
end
hold off ;
writematrix(CC,cat(2,nam,' arenas.csv'))
export_fig(cat(2,nam,' arenas.pdf'),'-painters')
close

bgim = imcomplement(bgim) ;
% Create output nam
%prm = sprintf('Width[%d,%d] minArea%d Length[%d,%d] cutoff%d', ...
%    minWidth, maxWidth, minArea, minLength, maxLength, cutoff) ; 

if ~isfolder(nam)
    mkdir (nam) ;    
end
outdir = cat(2,pwd,filesep,nam) ;
% Get arenas
CC = readmatrix(cat(2,nam,' arenas.csv')) ;
n = size(CC,1) ;    % no. of wells in the arena
track = cell(n,m) ; 
toc
tic    
if ~isfile(cat(2,nam,' X.csv')) % if output doesn't exists already  
    if isfile(cat(2,outdir,filesep,nam,sprintf(' fish%1.2d.mat',1)))
        fprintf('\nReading mat files from previous run\n')
        for j = 1:n
            load(cat(2,outdir,filesep,nam,sprintf(' fish%1.2d.mat',j)))
            ii = length(fish_track) ;                
            for i = 1:ii
                track{j,i} = fish_track(i) ;
            end
            fprintf('Read %d frames for fish %1.2d\n',ii,j) ;
        end
        img = read(V,m1+ii-1) ; % This is done to start from the right frame
        i = ii ; ik = 6 ;
        fprintf('\n%d frames (%d %%) left to analyze\n',m-ii,round(100*(m-ii)/m)) ;
    else
        img = read(V,m1-1) ; % This is done to start from the right frame
        i = 0 ; 
        if saveavi == 1
            ik = 1 ;
            v = VideoWriter(strcat(outdir,'\',nam,sprintf('-%1.2d',ik),'.avi')) ;    
            v.FrameRate = 2*V.FrameRate ;
            v.Quality = videoquality ;
            open(v); 
        end
    end
    while hasFrame(V) && i < m
        i = i+1;      
        img = mean(readFrame(V),3) ;
        img = img(150:1800,10:1070) ;
        if saveavi == 1 && ik <= 5
            IM = im2uint8(mat2gray(img)) ;
            rgb = cat(3, IM, IM, IM) ; %just for display
        end
        img = imcomplement(img) - bgim ; % Background subtraction
        img(img<0) = 0 ;
        img = im2uint8(mat2gray(img)) ; % image is rescaled between 0-255  
        imc = cell(n,1) ;
        for j = 1:n
            x0 = CC(j,4); x1 = x0+CC(j,6)-1 ;
            y0 = CC(j,3); y1 = y0+CC(j,5)-1 ;
            imc{j} = img(x0:x1, y0:y1) ; % sending from gpu to cpu
        end
        % analyze each well in the arena. Use parfor loop, but using too many
        % cores (>6) is not advantageous
        R = cell(n,1) ; C = cell(n,1) ;        
        parfor j = j1:j2 
            im = imc{j} ;
            idx = 0 ;
            [d1,d2] = size(im) ;
            prop = struct() ;
            prop.edges = nan ; 
            prop.img = nan ;
            prop.pos = nan ;
            prop.snr = nan ;
            prop.area = 0 ;
            % binarize image
            %thr = 255*graythresh(imc) ; %otsu's threshold
            
            if max(max(im)) > 100 % was 150
                thr = YenThresh(im)  ;           
                if border == 1
                    bwa = imclearborder(im > thr) ; 
                else
                    bwa = im > thr ;
                end                
                bwa = bwareafilt(bwa, [minArea,maxArea]) ;
                % I have tried breareaopen instead of bwareafilt, but
                % it leads to detection of multiple large edges
                
                % figure('Position',[1 1 800 800])
                % subplot(121)
                % imshow(im,[0 255])
                % subplot(122)
                % imshow(bwa, [])
                % title(sprintf('Frame %1.4d',i))
                % pause()
                % close()
                % while nnz(bwa) > maxArea
                %     thr = thr + 0.2*thr ; %increase threshold by 20%
                %     if border == 1
                %         bwa = imclearborder(imc{j} > thr) ; 
                %     else
                %         bwa = imc{j} > thr ;
                %     end
                %     bwa = bwareafilt(bwa, [minArea,maxArea]) ;
                % end
                % if nnz(bwa) < 2*minArea
                %     thr = thr - 0.2*thr ; %reduce threshold by 20%
                %     if border == 1
                %         bwa = imclearborder(imc{j} > thr) ; 
                %     else
                %         bwa = imc{j} > thr ;
                %     end
                %     bwa = bwareafilt(bwa, [minArea,maxArea]) ;
                % end
                
                %imwrite(im2uint8(bwa), sprintf('Step10 %s Well%d.tif', nam, j),'WriteMode','append')  
                bwa = bwpropfilt(bwa, 'MinorAxisLength', round([minWidth/px,maxWidth/px])) ; 
                bwa = bwpropfilt(bwa, 'MajorAxisLength', round([minLength/px,maxLength/px])) ;   
                %imwrite(bwa, sprintf('Step12 %s Frame%d Well%d.jpg', nam, i, j)) 
                bwa = imdilate(bwa,strel('diamond', 1)) ;

                % figure('Position',[1 1 800 800])
                % subplot(121)
                % imshow(im,[0 255])
                % subplot(122)
                % imshow(bwa, [])
                % title(sprintf('Frame %1.4d',i))
                % pause()
                % close
                
                S = regionprops(bwa,im,'Area','PixelIdxList','MeanIntensity','WeightedCentroid') ;
                obn = size(S,1) ; % no. of objects
                if  obn > 1  
                    iI = cat(1,S.MeanIntensity) ;
                    inp = floor(0.8*max(iI)) ;
                    i2 = find(iI>=inp) ;
                    
                    if ~isempty(i2)
                        ni2 = length(i2) ;
                    % for multiple bright objects, select the object
                    % with least ratio of perimeter to area                            
                        ss = regionprops(bwa,'Perimeter','Area') ;
                        s1 = cat(1,ss.Perimeter) ;
                        s2 = cat(1,ss.Area) ;                            
                        [ii,idx] = min(s1./s2) ;
                        
                        % If an object is close to edge, then dct will
                        % be small, but if it has high act, then the 
                        % product will tend to increase. On the otherhand,
                        % if an artifact is further from edge than the
                        % object, but its act (inverse aspect ratio) is lower, 
                        % then the product will tend to decrease
                        % [j3,i3] = max(dct*act) ;
                        % idx = i2(i3) ;
                    end
                    bwa = zeros(d1*d2,1) ;
                    bwa(S(idx).PixelIdxList) = 1 ;
                    bwa = reshape(bwa, [d1,d2]) ; 
                    S = regionprops(bwa,im,'Area','PixelIdxList','MeanIntensity','WeightedCentroid') ;                
                    
                    % figure('Position',[1 1 800 800])
                    % subplot(121)
                    % imshow(im,[0 255])
                    % subplot(122)
                    % imshow(bwa, [])
                    % title(sprintf('Frame %1.4d',i))
                    % pause()
                    % close
                end            
                
                if ~isempty(S) 
                    bg = im(~bwa) ; % Background pixels
                    SNR = S.MeanIntensity/std(double(bg))   ;   
                    ctd = S.WeightedCentroid ;     
                    prop.pos = px*ctd ; % mm
                    prop.snr = SNR ; 
                    prop.area = S.Area ;
                    if savemat == 1
                        imb = zeros(d1*d2,1) ;
                        im1 = reshape(imc{j},[d1*d2,1]) ;
                        imb(S.PixelIdxList) = im1(S.PixelIdxList) ;
                        imb = reshape(imb,[d1,d2]) ;
                        prop.img = sparse(imb) ;
                        if saveavi == 1 && ik <= 5
                            % Determine edges
                            bw = bwmorph(bwa,'remove') ;
                            bw = bwareafilt(bw,1) ;
                            edg_pts = struct2array(regionprops(bw, 'PixelList')) ;
                            R{j} = round(CC(j,4) + edg_pts(:,2) -1) ; 
                            C{j} = round(CC(j,3) + edg_pts(:,1) -1) ;
                            prop.edges = [R{j},C{j}] ;
                        end
                    end                 
                end            
            end
            track{j,i} = prop ;
        end
        if saveavi == 1 && ik <= 5
            for j = j1:j2
                for ij = 1:length(R{j})                        
                    rgb(R{j}(ij),C{j}(ij),1) = 255 ;
                    rgb(R{j}(ij),C{j}(ij),2) = 0 ;
                    rgb(R{j}(ij),C{j}(ij),3) = 80 ;
                end 
            end
            writeVideo(v, uint8(rgb)) ;
        end
        if mod(i,frame_buffer) == 0
            cs = round(i/toc,2) ;
            fprintf('%d of %d frames done for %s at %1.2f fps.', i, m, nam, cs)
            fprintf(' %1.2f minutes remain.\n',(m-i)/(cs*60))
            if saveavi == 1 && ik < 5
                close(v)                
                ik = ik+1 ;
                % Start new video object 
                v = VideoWriter(strcat(outdir,filesep,nam,sprintf('-%1.2d',ik),'.avi')) ;
                v.FrameRate = 2*V.FrameRate ; 
                v.Quality = videoquality;
                open(v);
            elseif saveavi == 1 && ik == 5
                close(v)                
                ik = ik+1 ;                    
            end  
        end
        if savemat == 1 && mod(i,10*frame_buffer) == 0
            fprintf('\n......Saving the mat files.......\n') ;
            for j = j1:j2
                fish_track = cat(1,track{j,:}) ;
                save(cat(2,outdir,filesep,nam,sprintf(' fish%1.2d.mat',j)), 'fish_track')    ;  
            end
        end
    end
    % if saveavi == 1 && mod(i,frame_buffer) ~= 0
    %     close(v)  % in case i is not a multiple of frame buffer
    % end   
    if savemat == 1 && mod(i,10*frame_buffer) ~= 0
        fprintf('\n......Saving the mat files.......\n') ;
        for j = j1:j2
            fish_track = cat(1,track{j,:}) ;
            save(cat(2,outdir,filesep,nam,sprintf(' fish%1.2d.mat',j)), 'fish_track')    ;  
        end
    end
    savetracks
end
toc
tic
% Generate heatmap by mapping the detected centroid
XY = zeros(max(CC(:,4))+(35/px),max(CC(:,3))+(35/px)) ;
cov = nan(n,1) ;
for j = 1:n                
    for i = 1:m
        if ~isnan(track{j,i}.pos)
            rX = ceil(track{j,i}.pos(1)/px) ; % mm -> pixel
            rY = ceil(track{j,i}.pos(2)/px) ; % mm -> pixel  
            if ~exist('pj','var')
                pj = double(full(track{j,i}.img)) ;
                [d1,d2] = size(pj) ;
            end
            % Transform well coordinates to plate coordinates for plotting            
            col = ceil(CC(j,3) + rX - 2) ;
            row = ceil(CC(j,4) + rY - 2) ;   
            XY(row,col) = XY(row,col) + 1 ;                
            pj = pj + double(full(track{j,i}.img)) ; 
        end
    end  
    if exist('pj','var')
        cov(j) = 100*nnz(pj>30)/(d1*d2) ;
    end
    if cov(j) < 1 %percent
        fprintf('\nWell %d is empty\n',j) ;   
    end
    clearvars pj ;
    toc
end
mx = max(max(XY))
th = prctile(XY(XY>1),[2 98])
imshow(XY,th) ; colormap(jetb); colorbar
export_fig(cat(2,nam,' heatmap.pdf'),'-painters','-transparent')
imwrite(uint16(XY),cat(2,nam,' heatmap.tif'))    
writematrix(cov,cat(2,nam,' coverage.csv'))
finish = datetime 