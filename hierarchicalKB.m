%% Author : Dr. Krishnashish Bose (Guo lab, UCSF)
%  Created on: 2024
%  Last modified on: 23-Jun-2025
%% This program reads .mat files containing neuronal activity matrices and
%% combines it with tailangle timeseries (read from a csv file)

clear; clc; close all;
f = filesep ;
fish = dir('*neurons.mat') ;
% fish = fish([fish.isdir])
dir0 = pwd
% fn = 22 ;
for fn = 1:length(fish)
    fname = erase(fish(fn).name,' neurons.mat')
    load(cat(2,fname,' neurons.mat')) ;
    n = length(neuron) ;        
    m = length(neuron{1}.Fflt) ;
    F = zeros(n,m) ;
    for i = 1:n
        F(i,:) = neuron{i}.Fflt ;
    end
    zF = zscore(F,[],2) ;
    F = zF - min(zF,2) ; % Makes all values positive
    
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
    
    cc = mat2gray(F) ;
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
end