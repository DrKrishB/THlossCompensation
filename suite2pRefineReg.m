%% Author: Dr. Krishnashish Bose (Guo Lab, UCSF)
%% Created on: 22-May-2020
%% Last modified on: 16-Jul-2024 
% clear;clc
% clearvars -except folder fn outdir
warning off
datetime

dir0 = pwd
f = filesep ;
pcname = strsplit(dir0,f) ;
pcname = pcname{2} ;
if isunix && strcmp(pcname,'wynton')
    nslots = getenv('NSLOTS') ;
    nslots = str2num(nslots) ;
    
    addpath(genpath('/wynton/home/guolab/krishb/scriptsKB')) ;
    movie_loc = '/wynton/scratch/krishb/' ; % temporary location with motion-corrected movies
else
    nslots = 1 ;
    fn = input('Which fish? ') ;
end
maxNumCompThreads(nslots) ;
if nslots > 2        
    parpool(nslots) ;
end
fish = dir('*dpf*') ;
fish = fish([fish.isdir])
fn
fname = fish(fn).name
if ~isfile (cat(2,fname,' neurons.mat'))
    prnt_detect = 1  % print overlay of detected neurons. Takes significantly more time.
    stepZ = 10  % step size in microns
    pixel = 0.977  % pixel size in microns
    min_dist = 3/pixel  % 5 um minimum allowed distance between neurons
    min_corr = 0.8  % for merging of wrongly splitted neurons)
    minArea = round(10/(pixel^2))  % minimum pixelcount of a neuron
    minLength = round(4/pixel)
    minWidth = round(3/pixel)
    maxEcn = 0.92  % maximum eccentricity of objects to qualify as cell

    options = CNMFSetParms()
    % options.temporal_parallel = false ;
    % options.method = 'spgl1'

    cd(fname)
    cd suite2p
    plane = dir('plane*') ; % Find all Z-planes
    nplane = length(plane)

    cnt = 0 ; % real neurons
    neuron_img = cell(nplane,1) ;
    projection = cell(nplane,1) ;
    activity_Braw = cell(nplane,1) ; % average mean intensity of each neuron (from pixels that I refined from suite2p's)
    activity_B = cell(nplane,1) ; % raw filtered Ca trace Fflt (Krishnashish)
    activity_Craw = cell(nplane,1) ; % raw fluorescence from suite2p (less noisy than my Fraw)
    activity_C = cell(nplane,1) ; % deconvolved Ca trace from CNMF/CaImAn
    activity_FS = cell(nplane,1) ; % Ca spikes from suite2p (very noisy and busy)
    activity_S = cell(nplane,1) ; % deconvolved Ca spikes from CNMF/CaImAn
    activity_snr = cell(nplane,1) ; % activity noise floor determined by CNMF/CaImAn

    for np = 1:nplane
        fld = strcat('plane',num2str(np-1)) ; % suite2p stores planes as Z-1
        cd (fld) ;
        if np <10
            outnam = cat(2, fname, ' ', 'Z0', num2str(np)) ;
        else
            outnam = cat(2, fname, ' ', 'Z', num2str(np)) ;
        end
        % Read suite2p output
        load('Fall.mat')
        n = length(stat) ;  % total neurons detected
        cd ..
        tic
        % Read motion-corrected movie
        mov = [] ;
        mvloc = strcat(movie_loc,'suite2pMCtif',f,fname,f,'suite2p',f,fld,f,'reg_tif',f,'*.tif') ;
        mvs = dir(mvloc) ;
        mcnt = length(mvs) ;
        if mcnt >= 1
	        mov = tiffreadVolume(cat(2,mvs(1).folder,f,mvs(1).name)) ;
        end
        if mcnt > 1
            for i = 2:mcnt
	            mov1 = mov ;
	            mov2 = tiffreadVolume(cat(2,mvs(i).folder,f,mvs(i).name)) ;
                mov = cat(3,mov1,mov2) ;
            end
        end
        [d1,d2,m] = size(mov)
        fprintf('\n%d motion-corrected movies of %s read in %d seconds\n',length(mvs),outnam,round(toc))
        tic
        if m == 1401
            mov(:,:,1) = [] ; % drop 1st frame
            m = size(mov,3)
        end
        % movmin = min(mov,[],3) ;  % Minimal projection
        % mov = mov - movmin ;
        % mov(mov<10) = 0 ;
        mov = mat2gray(mov) ; % entire movie normalized
        % I found that the movie size can be reduced to 0.45 GB from 1.36 GB by
        % saving mov as a sparse matrix with pixels above noise.
        % Compute projections
        movstd = mat2gray(imbilatfilt(std(mov,[],3))) ;  % Standard deviation projection
        movavg = imbilatfilt(mean(mov,3)) ;  % Maximal projection
        movmax = imbilatfilt(max(mov,[],3)) ;  % Maximal projection
        movmin = imbilatfilt(min(mov,[],3)) ;  % Minimal projection
        ref = movmax-movmin ;

        % Create a mask for entire fish
        IMF = medfilt2(imgaussfilt(ref,3),[10,10]) ;
        MASK = imbinarize(IMF,0.6*graythresh(IMF)) ;
        if length(find(MASK>0)) > 0.7*d1*d2
            MASK = imbinarize(IMF,0.9*graythresh(IMF)) ;
        end
        MASK = imdilate(MASK,strel('disk',2)) ;
        mps = find(MASK) ;

        %% Reject neurons based on morphology
        clear S
        thf = 0.1+0.01*np ;
        morph_filterKB    %Input stats | Output S
        n = length(S) ;
        fprintf('\n%d neurons detected in first round in %d minutes\n',n,round(toc/60)) ;
        tic
        %% Merge neurons based on overlap and similarity of activity
        clear P
        merge_neuronsKB   %Input S | Output P
        n = length(P) ;
        S1 = S ; S = P ;

        clear P
        merge_neuronsKB   %Input S | Output P
        n = length(P) ;
        S2 = S ;  S = P ;

        clear P
        merge_neuronsKB    %Input S | Output P
        n = length(P) ;

        %% Sort neurons based on position
        clear Q
        CC = zeros(k,2) ;
        for j = 1:n
            CC(j,:) = P{j}.center ;
        end
        [I,ord1] = sortrows(round(0.01*CC),[2 1]) ;
        q = I(:,2) ;
        qid = find(diff(q)==1) ;
        for ii = 1:length(qid)
            if mod(qid(ii),3) > 0
                q(qid(ii)) = q(qid(ii)) + 1 ;
            end
        end
        I(:,2) = q ;
        [I,ord2] = sortrows(I,[2 1]) ;
        order = ord1(ord2) ;
        CC = CC(order,:) ;
        for j = 1:n
            Q{j} = P{order(j)} ;
        end
        n = length(Q) ;
        clear R
        activity_filterKB ;
        n = length(R) ;

        fprintf('%d -> %d -> %d -> %d -> %d -> %d neurons refined in %d seconds\n', ...
            length(stat),length(S1),length(S2),length(S),length(Q),length(R),round(toc))
        tic
        % save temporal component as tif file
        B = zeros(n,m) ; FS = B ; craw = B ; Braw = B ;
        for j = 1:n
            B(j,:) = R{j}.Fflt ; % This is my own creation (Krishnashish)
            FS(j,:) = R{j}.Cflt ; % spikes from suite2p (very noisy and busy)
            Braw(j,:) = R{j}.Fraw ; % average mean intensity of each neuron (from pixels that I refined from suite2p's)
            craw(j,:) = R{j}.Craw ;  % raw fluorescence from suite2p (less noisy than my Fraw)
        end
        activity_B{np} = sparse(B) ; % filtered raw Ca trace Fflt (Krishnashish)
        activity_FS{np} = sparse(FS) ; % Ca spikes from suite2p (very noisy and busy)
        activity_Braw{np} = sparse(Braw) ; % average mean intensity of each neuron (from pixels that I refined from suite2p's)
        activity_Craw{np} = sparse(craw) ; % raw fluorescence from suite2p (less noisy than my Fraw)
        P = struct() ;
        P.p = 1 ;

        craw = zscore(craw,[],2) ; % for cnmfe's deconvolution
        craw = craw - min(craw,2) ; % no more negative values
        [C, parms, S] = deconv_temporal(craw, P, options) ;

        activity_C{np} = sparse(C) ; % deconvolved Ca trace from CNMF/CaImAn
        activity_S{np} = sparse(S) ; % deconvolved Ca spikes from CNMF/CaImAn
        activity_snr{np} = parms.neuron_sn ;

        nrn = zeros(d1,d2,n) ;
        for j = 1:n
            nrn(:,:,j) = full(R{j}.img) ;
            % imwrite(im2uint16(nrn(:,:,j)), cat(2, outnam, ' final neurons.tif'), 'WriteMode', 'append') ;
        end
        neuron_img{np} = sparse(reshape(nrn,[d1*d2,n])) ;

        % print the filtered detected neurons overlaid on Max-min projection
        nrnmax = reshape(mat2gray(max(nrn, [], 3)),[d1*d2,1]) ;
        Qnorm = quantilenorm(mat2gray(cat(2,ref,nrnmax,movstd)));
        projection{np} = Qnorm ;
        nam = cat(2, outnam, ' red_movmax green_nrnmax blue_movstd.tif') ;
        imwrite(im2uint16(reshape(Qnorm(:,1),[d1,d2])), nam) ;
        imwrite(im2uint16(reshape(Qnorm(:,2),[d1,d2])), nam, 'WriteMode','append') ;
        imwrite(im2uint16(reshape(Qnorm(:,3),[d1,d2])), nam, 'WriteMode','append') ;

        if prnt_detect == 1
        % print the filtered detected neurons overlaid on Max-Min projection
            if size(ref,2) == 1
                ref = reshape(ref, [d1,d2]) ;  % 2D Max-Min projection
            end
            max_val = 0.5*max(max(ref)) ;
            figure('Position', [1 1 800 500],'Visible','off');
            imshow(ref, [0 max_val]) ;            
            hold on;
            cn = zeros(n,2) ;
            for j = 1:n
                im = full(R{j}.img) ;
                bw = imbinarize(mat2gray(im),0.5) ;
                bw = bwareafilt(bw,1) ;
                cn(j,:) = R{j}.center ;
                %scatter(cn(j,1),cn(j,2), round(0.8*minArea), 'm', 'LineWidth',0.3) ;
                visboundaries(bw,'Color','m','LineWidth',0.2,'EnhanceVisibility',false) ;
                text(cn(j,1),cn(j,2),num2str(R{j}.ID),'Color','g','FontSize',1 ) ;
            end
            hold off
            export_fig (cat(2, outnam, ' neurons overlaid on Max-Min projection.pdf'), '-painters') ;
            close
        end
        % save(cat(2,outnam,' sorted.mat'), 'R')

        %% Append results from different Zplanes
        cnt = cnt +length(R) ;

        if np == 1
            neuron = R ;
            N = length(R) ;
        else
            neuron = [neuron, R] ; %concatenate structures
            N = cat(1,N,length(R)) ;
        end
        fprintf('%d neurons deconvolved in %s in %d minutes\n', n,outnam,round(toc/60)) ;
    end
    cd .. % get out of suite2p folder
    cd .. % get out of fish folder
    save(cat(2,fname,' neurons.mat'), 'neuron','neuron_img','projection','d1','d2','-v7.3')
    save(cat(2,fname,' neuron_activity.mat'), 'activity_Braw','activity_B','activity_Craw','activity_C','activity_FS','activity_S','activity_snr','-v7')
    % Save suite2p parameters
    writecell([{fname} N' sum(N)],cat(2,'NeuronCount.xlsx'),'WriteMode','append')
    parms = [{fname}, ops.fs, ops.nframes, ops.tau, ops.nimg_init, ops.batch_size, ...
            ops.smooth_sigma_time, ops.smooth_sigma, ops.snr_thresh, ...
            ops.maxregshiftNR, ops.pre_smooth, ops.sparse_mode, ...
            ops.diameter, ops.spatial_scale, ops.nbinned, ...
            ops.max_iterations, ops.threshold_scaling, ...
            ops.max_overlap, ops.allow_overlap, ops.spatscale_pix] ;
    if fn == 1
        parm_header = [{'movie'}, {'fs'}, {'nframes'}, {'tau'}, {'nimg_init'}, {'batch_size'}, ...
                {'smooth_sigma_time'}, {'smooth_sigma'}, {'snr_thresh'},...
                {'maxregshiftNR'}, {'pre_smooth'}, {'sparse_mode'}, ...
                {'diameter'}, {'diameter'}, {'spatial_scale'}, {'nbinned'}, ...
                {'max_iterations'}, {'threshold_scaling'}, ...
                {'max_overlap'}, {'allow_overlap'}, {'spatscale_pix'}] ;
        writecell(parm_header,'Suite2parms.xls')
    end
    writecell(parms,'Suite2parms.xls','WriteMode','append') ;    
else
    load(cat(2,fname,' neurons.mat')) ;
    n = length(neuron) ;
    m = length(neuron{1}.Fflt) ;
    F = zeros(n,m) ;
    for i = 1:n
        F(i,:) = neuron{i}.Fflt ;
    end
end
% z-score fluorescence
zF = zscore(F,[],2) ;
F = zF - min(zF,2) ; % Makes all values positive
cc = mat2gray(F) ;
sumF = sum(F,2) ; % summed across time

% Combine with tailangle data   
if isfile(cat(2,fname,' tailangle.csv'))
    taildata = readmatrix(cat(2,fname,' tailangle.csv')) ;
end
tailneuron = zeros(1,m) ;
Ts = linspace(taildata(1,1),taildata(end,1),m+1) ;
ii = find(taildata(:,1)<=Ts(2)) ;
tailneuron(1) = max(rmoutliers(abs(taildata(ii,2)),'percentiles',[0 90])) ;
for i = 2:m
    ii = find(taildata(:,1)>Ts(i) & taildata(:,1)<=Ts(i+1)) ;
    tailneuron(i) = max(rmoutliers(abs(taildata(ii,2)),'percentiles',[0 90])) ;
end
Ts = linspace(taildata(1,1),taildata(end,1),m) ;
save(cat(2,fname,' tailneuron.mat'),'tailneuron','Ts','F','-v7')

tailneuron(tailneuron<5) = 0 ; % angles less than 5 degree are ignored 
cc(n+1,:) = mat2gray(tailneuron) ; % Add tail data as an extra row to neural data

% Perform hierarchical clustering
tic
D = pdist(cc, 'spearman') ; % I tried GPU. It's not much faster.
tree = linkage(D, 'average') ;    
fprintf('\nHierarchical clustering completed in %d minutes\n',round(toc/60)) ;
nc = n ;
t = 2 ; % starting no. of clusters  
while t < 300  % the no. of neurons belonging to tail cluster should not suddenly drop by more than 20% of total neurons  
    t = t+1 ;
    clust = cluster(tree,'maxclust',t) ;
    for i = 1:t
        idx = find(clust==i) ;            
        if ismember(n+1,idx)
            fprintf('%d neurons belonging to cluster %d (of %d) are correlated to tail movements\n',length(idx),i,t) ;
            nc = [nc length(idx)] ;
            break ;
        end
    end 
    if length(nc) > 3 && length(idx) < 0.8*n
        % average difference in the last 10 nc is less than 10 neurons
        nc = nc(2:4) ; % last 3 tailneuron count
        if max(abs(diff(nc))) < 10
            fprintf('\nOptimal clusters reached')
            t = t-2 % optimal cluster reached 
            break ;
        end
    end         
end  
% go back if too few neurons
while nc(1) < 0.4*n
    t = t-1 
    clust = cluster(tree,'maxclust',t) ;
    for i = 1:t
        idx = find(clust==i) ;            
        if ismember(n+1,idx)
            fprintf('%d neurons belonging to cluster %d (of %d) are correlated to tail movements\n',length(idx),i,t) ;
            nc = length(idx) ;
            break ;
        end
    end 
end
clust = cluster(tree,'maxclust',t) ;
clusters = cell(t,1) ;
nclust = zeros(t,1) ;
for i = 1:t
    clusters{i} = find(clust==i) ;
    nclust(i) = length(clusters{i}) ;
    if ismember(n+1,clusters{i})
        fprintf('\nSelected %d neurons belonging to cluster %d (of %d) correlated to tail movements\n',nclust(i),i,t) ;
        tailneuron_cluster = i ;
    end
end  
save(cat(2,fname,' tailneuron.mat'),'tailneuron_cluster','clusters','nclust','-append')   
datetime
%% Registration of neurons 
% Additional inputs: Z-stack green and red
pos = zeros(n,3) ;
for i = 1:n
    pos(i,:) = round(neuron{i}.xyz) ;
    pos(i,3) = pos(i,3) + stepZ/2 ; % for continuous Z-ramping
end     
if isunix && strcmp(pcname,'wynton')
    pd = '/wynton/home/guolab/krishb/registration/' ; % Primary Directory
    ref = strcat(pd,'reg2 211016 7dpf HRmtzl03 green rotate0.8 [0.7 0 0].tif') ;
elseif ~isunix
    % ref = ['C:\Users\Krishnashish\OneDrive\Programs\scriptsKB\registration\reg2 211016 7dpf HR' ...
    ref = ['D:\OneDrive\scriptsKB\registration\reg2 211016 7dpf HRmtzl03 green rotate0.8 [0.7 0 0]' ...
        '\reg2 211016 7dpf HRmtzl03 green rotate0.8 [0.7 0 0].tif'] ;
    % pyenv('Version', 'C:\Users\Krishnashish\anaconda3\python.exe') ;
end

%% ANTs registration with Zbrain
if ~isfile(cat(2,fname,' reg2zbrain_XYZ.csv'))
    np = py.importlib.import_module('numpy') ;
    ants = py.importlib.import_module('ants') ;
    pos = np.array(pos) ;
    pts = pos.tolist() ;
    
    % Find warping transformation
    ref = ants.image_read(ref) ;
    ants.set_spacing(ref, {0.65, 0.65, 2.0}) ;
    sz = double(ref.shape) ;
    dY1 = sz(2) ;   dX1 = sz(1) ;   dZ1 = sz(3) ;    
    [dY1,dX1,dZ1]
    
    if ~isfile(strcat(dir0,f,fname,' reg2affine.mat'))
        tic
        stack = ants.image_read(strcat(fname," stack green mask.tif")) ;
        ants.set_spacing(stack, {0.977, 0.977, 2.0}) ;
    
        reg = ants.registration( fixed=ref , moving=stack, type_of_transform='SyNAggro') ;        
        fdTr = cellstr(string(reg{'fwdtransforms'}))     
    
        % save registered green channel
        trxG = reg{'warpedmovout'} ;
        trxG = trxG.numpy() ;
        trxG = double(trxG) ;
        stack = zeros(dY1,dX1,dZ1) ;
        for ij = 1:dZ1
            stack(:,:,ij) = trxG(:,:,ij)' ;
        end
        stack = uint16(stack) ;
        imwrite(stack(:,:,1),strcat(fname,' reg2zbrain_green.tif')) ;
        for ij = 2:dZ1
            imwrite(stack(:,:,ij),strcat(fname,' reg2zbrain_green.tif'),'WriteMode','append') ;
        end
        
        % Align red channel and save
        red = ants.image_read(strcat(fname," stack red.tif")) ;
        ants.set_spacing(red, {0.977, 0.977, 2.0}) ;
        trxR = ants.apply_transforms(fixed = ref, moving = red, transformlist = fdTr) ;
        trxR = trxR.numpy() ;
        trxR = double(trxR) ;
        red = zeros(dY1,dX1,dZ1) ;
        for ij = 1:dZ1
            red(:,:,ij) = trxR(:,:,ij)' ;
        end
        red = uint16(red) ;
        imwrite(red(:,:,1),strcat(fname,' reg2zbrain_red.tif')) ;
        for ij = 2:dZ1
            imwrite(red(:,:,ij),strcat(fname,' reg2zbrain_red.tif'),'WriteMode','append') ;
        end
        fprintf('\nRegistration completed in %d minutes using %s and %s\n',round(toc/60),fdTr{1},fdTr{2}) ;
        
        movefile(fdTr{1},strcat(dir0,f,fname,' reg2fwdTr.gz')) ;
        movefile(fdTr{2},strcat(dir0,f,fname,' reg2affine.mat')) ;
    end
    fdTr{1} = strcat(dir0,f,fname,' reg2fwdTr.gz') ;
    fdTr{2} = strcat(dir0,f,fname,' reg2affine.mat') ;

    % Align points stored in csv file
    XYZ = zeros(n,3) ;
    tx = ants.read_transform(fdTr{2}) ;
    rx = tx.invert() ;          
    for i = 1:n
        pt = rx.apply_to_point(pts{i}) ;
        pt = cell2mat(cell(pt)) ;
        XYZ(i,:) = double(pt) ;
    end
    XYZ = round(XYZ) ;
    pos = round(double(pos)) ;
    writematrix([pos XYZ],cat(2,fname,' reg2zbrain_XYZ.csv')) ;  
else
    XYZ = readmatrix(cat(2,fname,' reg2zbrain_XYZ.csv')) ;
    XYZ = XYZ(:,4:6) ;
end

%% Save registered points as neuron3D matrix
if ~isfile(cat(2,fname,' neurons3Dreg.mat'))
    tic
    Yd = 500 ; %um
    Xd = 1000 ; %um
    Zd = 200 ; %x2um
    neuron3D = cell(n,1) ; %sparse 3D stack for each neuron
    x = round(XYZ(:,1)) ;
    y = round(XYZ(:,2)) ;
    z = round(0.5*XYZ(:,3)) ; % 1 um -> 2 um z-spacing
    parfor i = 1:n
        nrn3d = false(Yd,Xd,Zd) ;  
        % Each neuron has 5x5x3 =  75 voxels
        for i1 = -2:2 % 1 um -> 5 um
            xx = min(Xd,max(1,x(i)+i1)) ;
            for i2 = -2:2 % 1 um -> 5 um
                yy = min(Yd,max(1,y(i)+i2)) ;
                for i3 = -1:1 % 2 um -> 6 um
                    zz = min(Zd,max(1,z(i)+i3)) ;
                    nrn3d(yy,xx,zz) = true ;
                end
            end
        end
        neuron3D{i} = sparse(reshape(nrn3d,[Yd*Xd*Zd,1])) ;
    end
    fprintf('\n %d neuron3D cells created in %d minutes\n',n,round(toc/60))   
    
    save(cat(2,fname,' neurons3Dreg.mat'),'neuron3D','sumF','-v7')
else
    load(cat(2,fname,' neurons3Dreg.mat'))
end
if ~exist('ovpmask','var')
    % Mask assignment
    tic
    if isunix && strcmp(pcname,'wynton')
        load('/wynton/home/guolab/krishb/registration/CombinedMasks701reg2.mat') ; % Anatomical masks
    else
        load('CombinedMasks701reg2.mat') % I would have already added to path
    end
    nmask = length(maskname) ;
    mid = cell(n,1) ; % assigned masks
    ovpmask = zeros(nmask,n) ; % overlapping pixels
    parfor j = 1:nmask
        for i = 1:n
            ovpmask(j,i) = nnz(neuron3D{i} & bwmask{j}) ; % apply AND operation to the two 3D logical matrices        
        end
    end
    ovpmask = sparse(ovpmask) ;
    fprintf('\n %d x %d overlap matrix computed in %d minutes\n',nmask,n,round(toc/60))  
    save(cat(2,fname,' neurons3Dreg.mat'),'ovpmask','-append')
else
    nmask = size(ovpmask,1) ;
end
bn = ovpmask >= 8 ; % At least 10% voxels should overlap with the mask
nn = length(find(sum(bn)==0)) ; % neurons not assigned to any mask
fprintf('\n%d (%1.1f%%) of %d neurons remain unlabeled in %s\n',nn,round(100*nn/n,1),n,fname) ;

%% Find neuron count and activity in different anatomical regions   
% writecell({'NeuronID','X_um','Y_um','Z_um','sumF','ZscoredF','sumOverlap','sumLabels','Labels'}, cat(2,fname,' neurons.xlsx'))
% writematrix(cat(2,(1:n)',x,y,2*z,round(F,3),round(zF,3),sum(np)',sum(bn)'),cat(2,fname,' neurons.xlsx'),'Range','A2')
A = zeros(1,nmask) ;    % Activity of each mask
N = zeros(1,nmask) ;    % Neurons in each mask
for i = 1:n % neurons
    J = find(bn(:,i)) ; % Assigned masks
    for k = 1:length(J)
        A(J(k)) = A(J(k)) + sumF(i) ; % Add neuron fluorescence to mask
        N(J(k)) = N(J(k)) + 1 ; % Add neuron fluorescence to mask
    end
end
if max(A) < 100
    A = round(A,3) ;
else
    A = round(A) ;
end
writematrix(A,'Activity per mask reg2zbrain_old.xlsx','Range',sprintf('A%d',fn)) ;
writematrix(N,'Neurons per mask reg2zbrain_old.xlsx','Range',sprintf('A%d',fn)) ;



