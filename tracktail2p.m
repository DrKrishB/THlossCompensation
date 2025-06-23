%% Author : Dr. Krishnashish Bose (Guo lab, UCSF)
%  Created on: 01-Oct-2020
%  Last modified on: 23-Jun-2025

%% The working directory needs to be a folder containing the tail tracking avi fishs.

clear; clc; close all;
warning off ;

bd = 5 ;  %extra padding on the left-right of base & tip points
% Tend = 1445 ; %Last timepoint obtained from ThorImageLS
% tailinfo = 'E:\TH Ca++ Tailmovies\tail info.xlsx' ;
% fn = input('File no. : ') ; % in the excel sheet, not folder
% 
% fname = char(info.FishName(fn))
% start = info.Start(fn) ; 
% stop = info.Stop(fn) ;
% fps = info.FrameRate(fn) ;
% m = stop - start + 1 ;
dir0 = pwd ;
fish = dir('*.avi') ;  %get fishname 
for i = 1:length(fish)
    disp(fish(i).name);
end
if ~exist('fn','var')
    fn = input('\nEnter fish index to analyze : ','s') ;
end
fn = str2num(fn) ;
j = fn(1) ;
str = strsplit(erase(fish(j).name,'.avi')) ;
str(end) = [] ;
fname = strjoin(str)

avi = VideoReader(fish(j).name) 
d1 = avi.Width; d2 = avi.Height ;
m = avi.NumFrames ; fps = avi.FrameRate ;
tailbw = cell(m,1) ;
Tend = 0 ;
if ~isfile("tailinfo.xlsx")
    writecell({'Filename','NumFrames','FrameRate','baseX','tipX','baseY', ...
        'tipY','Rotate','TailRest','Thresh','Start2p','End2p'},'tailinfo.xlsx')
end
info = readtable("tailinfo.xlsx") ;

if isempty(info.baseX)
    detect_tail
elseif length(info.baseX) < j
    detect_tail
else
    base_point(1) = info.baseX(j) ;
    tip_point(1) =  info.tipX(j) ;
    base_point(2) = info.baseY(j) ;
    tip_point(2) = info.tipY(j) ;
    R = info.Rotate(j) ;
    tail_rest = info.TailRest(j) ;
    th = info.Thresh(j) ;
    tail_length = round(((base_point(1) - tip_point(1))^2 + (base_point(2) - tip_point(2))^2)^0.5) ;
end

i = 1 ; tic
if ~isfile(cat(2,fname,' tailskeleton.mat'))
    for j = 1:length(fn)
        avi = VideoReader(fish(fn(j)).name) ;
        % if exist('start','var')
        %     img = read(avi,start-1) ;
        % else
        %     start = 1 ;
        % end    
        while hasFrame(avi) 
            if R ~= 0
                img = imrotate(mean(readFrame(avi),3), R) ;
            else
                img = mean(readFrame(avi),3) ;
            end
            imf = imgaussfilt(img,2)-imgaussfilt(img,10) ;
            imc = imf(:,base_point(1)-bd:tip_point(1)+bd) ; 
            %mp = max(max(imc)) ;
            
            bw = imfill(imbinarize(imc, th),'holes') ;
            bw = bwpropfilt(bw,'MajorAxisLength', [ceil(0.3*tail_length), Inf]) ;
            L = labelmatrix(bwconncomp(bw)) ;
            S = struct2array(regionprops(bw, 'MajorAxisLength')) ;
            [I,J] = max(S) ;
            bw = ismember(L,J) ; 
            ar = struct2array(regionprops(bw,'FilledArea')) ;
            if ar > tail_length*30
                bw = imfill(imbinarize(imc, 1.1*th),'holes') ;
                bw = bwpropfilt(bw,'MajorAxisLength', [ceil(0.3*tail_length), Inf]) ;
                L = labelmatrix(bwconncomp(bw)) ;
                S = struct2array(regionprops(bw, 'MajorAxisLength')) ;
                [I,J] = max(S) ;
                bw = ismember(L,J) ; 
            end   
            if gpuDeviceCount > 0
                sk = gather(bwmorph(gpuArray(bw),'thin',Inf)) ; % fast, but need GPU 
            else
                sk = bwskel(bw,'MinBranchLength', ceil(0.2*tail_length)) ; % slow
            end
            tailbw{i} = sparse(sk) ;    
            if mod(i,1000) == 0
                fprintf('%d seconds done at %1.2f fps\n', round(i/fps), round(i/toc,2))        
            end
            if mod(i,100000) == 0
                save(cat(2,fname,' tailskeleton.mat'),'tailbw','-v7.3') ;
            end
            i = i+1 ;
        end        
    end
    Tend = Tend + avi.Duration     
    toc 
    save(cat(2,fname,' tailskeleton.mat'),'tailbw','-v7.3') ;
else
    load(cat(2,fname,' tailskeleton.mat')) ;
end
    
m = length(tailbw) % This may be bigger than initial allocation
T = linspace(0,Tend,m)' ;
ns = 8 ; %Divide tail into 'ns' equal segments
%points representing each segment
xs = cell(m,ns) ;   
ys = cell(m,ns) ;
%slope of each segment 
tailangle = nan(m,ns) ; 
% Frame size
[d1,d2] = size(tailbw{1}) ;
if ~isfile(cat(2,fname,' tailpoints.mat'))
    tic
    for i = 1:m
        bw = full(tailbw{i}) ;  
        pts = cell2mat(bwboundaries(bw)) ;  %the points are sorted along the contour 
        %the 1st column is y & 2nd column is x
        if length(pts) > ns
            np = round(0.5*length(pts)) ;
            pts = pts(1:np,:) ;   
            %find the correct starting point
            %whichever end-point is closer to centre
            x0 = 0; y0 = round(d1/2) ;
            x1 = pts(1,2) ; y1 = pts(1,1) ;
            dst1 = (x0-x1)^2 + (y0-y1)^2 ;
            x1 = pts(np,2) ; y1 = pts(np,1) ;
            dst2 = (x0-x1)^2 + (y0-y1)^2 ;
            if dst2 < dst1
                pts = flip(pts) ;
            end
            x = pts(:,2) ; y = pts(:,1) ;
            L = round(np/ns) ; %length of each segment  
            i0 = 1 ;
            for ij = 1:ns            
                if ij < ns
                    i1 = i0+L-1 ;
                else
                    i1 = np ;
                end
                xp = x(i0:i1) ;
                yp = y(i0:i1) ;
                xr = max(xp) - min(xp) ;
                yr = max(yp) - min(yp) ;
                if xr >= yr
                    pp = polyfit(xp,yp,3) ; %3rd order fit 
                    xs{i,ij} = xp ;
                    ys{i,ij} = polyval(pp,xp) ;
                else
                    pp = polyfit(yp,xp,3) ; %3rd order fit 
                    ys{i,ij} = yp ;
                    xs{i,ij} = polyval(pp,yp) ;
                end 
                %get slope by 1st order fitting
                ps = polyfit(xs{i,ij},ys{i,ij},1) ; % 1st order fit  to get slope  
                tailangle(i,ij) = -round(atand(ps(1)),1) ;% 1st coefficient is slope
    
                i0 = i0+L ;
            end
        end
        if mod(i,100000) == 0
            fprintf('%d seconds tracked at %1.2f fps\n', round(i/fps), round(i/toc,2))      
            writematrix(cat(2, T, tailangle),cat(2,fname,' tailtrack.csv'))
        end
    end
    writematrix(cat(2, T, tailangle),cat(2,fname,' tailtrack.csv'))
    save(cat(2,fname,' tailpoints.mat'), 'xs','ys')
    toc
else
    load(cat(2,fname,' tailpoints.mat')) ;
    data = readmatrix(cat(2,fname,' tailtrack.csv')) ;
    T = data(:,1) ;
    tailangle = data(:,2) ;
end

%% Get one angle
x0 = base_point(1)-bd ;
y0 = round(d1/2) ;
theta = zeros(m,1) ;
for i = 1:m
    xps = x0+cat(1,xs{i,4:7}) ;
    yps = cat(1,ys{i,4:7}) ;
    xc = mean(xps) ; yc = mean(yps) ;
    theta(i) = -atand((yc-y0)/(xc-x0)) ; 
end
theta = fillmissing(theta,'nearest') ; %remove any nan values
% remove baseline in a 1-minute window
tt = smoothdata(theta,'movmedian',round(60*fps)) ; %recording duration was 24 minutes
theta = theta-tt ;    
writematrix(cat(2,T,theta),cat(2,fname,' tailangle.csv')) ;
figure('Position', [ 10 10 1900 200 ])
plot(T,theta,'k-') ; axis tight ;
ylim([-90,90])
xlim([0,Tend])
xlabel('Time (seconds)')
ylabel('Tail angle (degrees)')
xticks(0:30:1450)
yticks(-90:30:90)
set(gca,'TickLength',[.001 .002])
title(replace(fname,'_',' '))
export_fig (cat(2,fname,' tailangle.pdf'),'-painters','-transparent')    
close
toc 

%% Print projections
tic
winT = 5; %seconds
cl = [[1 0 1];[1 0.6 0];[0 1 1];[0.7 0.2 1];[0 1 0];[1 0 0];[0 0 1];[0.6 0 0.2]] ;
win = ceil(winT*fps) ; %No. of frames in each projection of 'winT' seconds	
%Create figure
rows = 7 ; cols = 14 ;
figure('Position', [1 1 150*cols 150*rows]) 	
ha = tight_subplot(rows,cols, [.01 .01],[.01 .01],[.01 .01]) ;  
i = 1 ; %Frame counter
ik = 1; %Plot counter
k = 0 ; %PDF page counter 
for jj = 1:length(fn)
    avi = VideoReader(fish(fn(jj)).name) ;    
    mm = avi.NumFrames ;
    on = 1 ;
    while hasFrame(avi) && i < m
        off = min(mm,on+win-1) ;
        % mov = squeeze(mean(read(avi,[start+on-1,start+off-1]),3)) ;    
        mov = squeeze(mean(read(avi,[on,off]),3)) ;    
        axes(ha(ik)); 
        img = imrotate(max(mov,[],3),R) ;
        img = img - imgaussfilt(img,10) ;
        img(img<0) = 0 ;
        %imwrite(uint8(img), cat(2,fname,' projections.tif'), 'WriteMode', 'append')
        imshow(mat2gray(img), [0 0.5]) ; 
        hold on ; 
        xrange = zeros(ns,1) ;
        yrange = zeros(ns,1) ;
        for j = 1:ns
            xps = x0 + cat(1,xs{i:min(m,i+win-1),j}) ;
            if ~isempty(xps)
                yps = cat(1,ys{i:min(m,i+win-1),j}) ;
                scatter(xps,yps,1,cl(j,:)) ;
                xrange(j) = max(xps) - min(xps) ;
                yrange(j) = max(yps) - min(yps) ;
            end
        end
        hold off;   
        xL = round(mean(xrange(1:4))) ;
        xR = round(mean(xrange(5:ns))) ;
        yL = round(mean(yrange(1:4))) ;
        yR = round(mean(yrange(5:ns))) ;
        title(sprintf('%d - %d s | %d,%d,%d,%d', round(T(i)),round(T(min(m,i+win-1))),xL,xR,yL,yR),'FontSize',7)                
        if mod(ik,rows*cols) == 0
            ik = 1 ; k = k+1;
            export_fig (cat(2,fname, ' projections',num2str(k),'.pdf'), '-painters', '-transparent') ;
            close
            figure('Position', [1 1 150*cols 150*rows]) 
            ha = tight_subplot(rows,cols, [.01 .01],[.01 .01],[.01 .01]) ; 
        else
            ik = ik+1;
        end
	    i = i+win ;
        on = on+win ;
    end
end
if ik > 1
    k = k+1;
    export_fig (cat(2,fname, ' projections',num2str(k),'.pdf'), '-painters', '-transparent') ;
    close 
end 
toc
exit