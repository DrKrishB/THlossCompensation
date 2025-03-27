%% You just input the mask. The program reads *neuron3d.mat data from both 
% control and defective groups of larvae and outputs the p-value and 
% percentage difference in neuron activity, count and firing rate.
clear; clc
%function [pA,ppA,mA,pC,ppC,mC] = hyperactivity_mapZebrain(maskloc)
ot = 10 ; % minimum overlap between neuron and mask
dir1 = 'X:\suite2p_filtered_2025Mar\' ;
dY = 500; dX = 1000; dZ = 200 ;
% 1st mask
maskloc = 'C:\Users\Krishnashish\OneDrive\scriptsKB\registration\HRmtzl03rotds\mapZebrain_stacksKBreg\' ;
% maskname = 'GAD1b_GFP' ;
% maskname = 'vglut2a_dsRed' ;
% maskname = 'drd2aGal4' ;
maskname = 'chataGal4GFP' ;
% maskname = 'reticulospinal_backfills' ;

tif = mat2gray(tiffreadVolume(cat(2,maskloc,maskname,'.tif'))) ;
mask = tif > 0.5 ; % convert to logical
sz = size(mask) 
if nnz(mask) > 0.01*sz(1)*sz(2)*sz(3)
    mask = tif > 0.7 ;
elseif nnz(mask) < 1000
    mask = tif > 0.3 ;
end
if sz(1) == dY && sz(2) == dX && sz(3) == dZ
    mask = sparse(reshape(mask,[sz(1)*sz(2)*sz(3),1])) ;
else % mask need to be resized
    mask = imresize3(mask,[dY,dX,dZ],"Antialiasing",false) ;
    im = im2uint16(mask) ;
    imwrite(im(:,:,1),sprintf('Resized %s.tif',maskname))
    for z = 2:200
        imwrite(im(:,:,z),sprintf('Resized %s.tif',maskname),'WriteMode','append') ;
    end
    mask = sparse(reshape(mask,[dY*dX*dZ,1])) ;
end

% 2nd mask
maskloc2 = 'C:\Users\Krishnashish\OneDrive\scriptsKB\registration\CombinedMasks701reg2new\' ;
maskname2 = 'MLR_KB' ;
% maskname2 = '701 glut_subpallium' ;

tif2 = mat2gray(tiffreadVolume(cat(2,maskloc2,maskname2,'.tif'))) ;
mask2 = tif2 > 0.5 ; % convert to logical
sz = size(mask2) 
if nnz(mask2) > 0.01*sz(1)*sz(2)*sz(3)
    mask2 = tif > 0.7 ;
elseif nnz(mask2) < 1000
    mask2 = tif > 0.3 ;
end
if sz(1) == dY && sz(2) == dX && sz(3) == dZ
    mask2 = sparse(reshape(mask2,[sz(1)*sz(2)*sz(3),1])) ;
else % mask2 need to be resized
    mask2 = imresize3(mask2,[dY,dX,dZ],"Antialiasing",false) ;
    im = im2uint16(mask2) ;
    imwrite(im(:,:,1),sprintf('Resized %s.tif',mask2name))
    for z = 2:200
        imwrite(im(:,:,z),sprintf('Resized %s.tif',mask2name),'WriteMode','append') ;
    end
    mask2 = sparse(reshape(mask2,[dY*dX*dZ,1])) ;
end

% combined mask
mask = mask & mask2 ;
maskname = strcat(maskname,'&',maskname2)
mvol = nnz(mask) 
if mvol > ot
    dir0 = fileparts(pwd) ;    
    fishA = dir('*ctrl*neurons3Dreg.mat') ;
    nA = length(fishA) 
    A = zeros(1,nA) ;
    cntA = zeros(1,nA) ;
    spkA = zeros(1,nA) ;
    tic
    for fn = 1:nA    
        neuron3D = [] ; sumF = [] ;
        load(fishA(fn).name) ; 
        fname = erase(fishA(fn).name,' neurons3Dreg.mat') ;
        load(cat(2,dir1,fname,' neuron_activity.mat'),'activity_S')
        s = cat(1,activity_S{:}) ;
        %tc = tailneuron>5 ; 
    
        n = length(neuron3D) ; % no. of neurons
        S = zeros(n,1) ;
        nv = zeros(n,1) ; % no. of overlapping voxels
        parfor i = 1:n
            S(i) = length(findpeaks(s(i,:))) ;
            nv(i) = nnz(mask & neuron3D{i}) ;
        end    
        J = find(nv >= ot) ; % At least 'ot' voxels should overlap with the mask
       
        cntA(fn) = length(J) ;
        A(fn) = round(sum(sumF(J)),2) ;  
        spkA(fn) = round(sum(S(J))/length(J),2) ;
    end
    toc
    fishB = dir('*mtzh*neurons3Dreg.mat') ;
    nB = length(fishB) 
    B = zeros(1,nB) ;
    cntB = zeros(1,nB) ;
    spkB = zeros(1,nB) ;
    tic
    for fn = 1:nB
        neuron3D = [] ; sumF = [] ;
        load(fishB(fn).name) ;
        fname = erase(fishB(fn).name,' neurons3Dreg.mat') ;
        load(cat(2,dir1,fname,' neuron_activity.mat'),'activity_S')
        s = cat(1,activity_S{:}) ;
        %tc = tailneuron>5 ;
        
        n = length(neuron3D) ; % no. of neurons
        S = zeros(n,1) ;
        nv = zeros(n,1) ; % no. of overlapping voxels
        parfor i = 1:n
            S(i) = length(findpeaks(s(i,:))) ;
            nv(i) = nnz(mask & neuron3D{i}) ;
        end    
        J = find(nv >= ot) ; % At least 'ot' voxels should overlap with the mask
         
        cntB(fn) = length(J) ;
        B(fn) = round(sum(sumF(J)),2) ; 
        spkB(fn) = round(sum(S(J))/length(J),2) ;
    end
    toc
    
    [pA,mA] = MannWhitney(A,B) 
    ppA = permtest(A,B) 
    [pC,mC] = MannWhitney(cntA,cntB) ;
    ppC = permtest(cntA,cntB) ;
    [pS,mS] = MannWhitney(spkA,spkB) ;
    ppS = permtest(spkA,spkB) ;
    
    savnam = sprintf('%s Hyperactivity bn%d %s.xlsx',maskname,ot,datestr(now,'yymmdd_HHMM'))
    writematrix([pA ppA mA nan A nan B], savnam)
    writematrix([pC ppC mC nan cntA nan cntB], savnam,'Sheet',2)
    writematrix([pS ppS mS nan spkA nan spkB], savnam,'Sheet',3)
end