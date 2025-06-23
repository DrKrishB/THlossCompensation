%% Author : Dr. Krishnashish Bose (Guo lab, UCSF)
%  Created on: Mar-2025
%  Last modified on: 23-Jun-2025

% Let's print clusters. Navigate to folder containing *neurons3Dreg.mat,
% which contain 'ovpmask', a 2D matrix of neuronal overlap with masks.
% To incorporate newer masks, use 'update_ovpmask.m'
% imx is the upper limit of raster plot 
% activity of tail-neuron cluster (red) and non-tail-neurons (blue) are
% divided by the total no. of neurons in that mask
% pmx is the upper limit of activity plot
% Usage: tailNeuron2x('all','j1',j) or tailNeuron2x('mtzh','j1',j)
% If 'j2' is specified, neurons common to both j1 and j2 are printed
function [] = tailNeuron2x(varargin)
dir0 = fileparts(pwd) ;
warning off ; tic
pcname = strsplit(dir0,filesep) ;
pcname = pcname{2} ;
if isunix && strcmp(pcname,'wynton')    
    addpath(genpath('/wynton/home/guolab/krishb/scriptsKB')) ;
end
if contains(pwd,'mapZebrain')
    load('mapZebrain.mat','maskname')
else
    % load('CombinedMasks707reg2.mat','maskname')
    load('ExpandedMasks.mat','maskname')
end

load('parulaRed.mat') ; % custom colormap by KB

imx = 10 ; % maximum colormap scale for raster plots
ot = 10 ; % minimum overlap between neuron and mask

% default parameter-values
prnt = 1 ;  
if ~isunix
    ext = 'pdf' ;
else
    ext = 'eps' ; 
end
pmx = 2 ;
switch varargin{1}  
    case 'ctrl'
        fish = dir('*ctrl*neurons3Dreg.mat') ;
        outdir = sprintf('Control TailNeuronClusters_%s',datestr(now,'yymmdd')) ;
    case 'mtzh'
        fish = dir('*mtzh*neurons3Dreg.mat') ;
        outdir = sprintf('Defective TailNeuronClusters_%s',datestr(now,'yymmdd')) ;
    case 'all'
        fish = dir('*neurons3Dreg.mat') ;
        outdir = sprintf('TailNeuronClusters_%s',datestr(now,'yymmdd')) ;
end
fn1 = 1 ;
fn2 = length(fish) ;

% Process name-value arguments to overwrite default values
for i = 2:2:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}  
            case 'j1'
                j1 = varargin{i+1} ;
                j2 = j1 ;
            case 'j2'
                j2 = varargin{i+1} ;
            case 'fn1'
                fn1 = varargin{i+1} ;
                fn2 = fn1 ;
            case 'fn2'
                fn2 = varargin{i+1} ;                
            case 'pmx' % upper limit of total activity plot
                pmx = varargin{i+1} ;
            case 'pdf'
                ext = 'pdf' ;            
        end
    end
end
mkdir(outdir)
outdir = cat(2,pwd,filesep,outdir,filesep) ;
if j1 ~= j2
    mnam = cat(2,maskname{j1},' & ',maskname{j2})
else
    mnam = maskname{j1} ;
end
mnam = replace(mnam,'Diencephalon ','Die ') ;
mnam = replace(mnam,'Mesencephalon ','Mes ') ;
mnam = replace(mnam,'Rhombencephalon ','Rho ') ;
mnam = replace(mnam,'Telencephalon ','Tel ') ;

mnam = cat(2,num2str(j1),' ',mnam) ;
savnam = cat(2,outdir,mnam,'.xlsx') 


writecell({'Fishname','Neurons','Tailneurons','TailFraction', ...
    'TaildfF','TailCorr','NonTaildfF','NonTailCorr','TotalClusters', ...
    'Cluster01','Cluster02', 'Cluster03', 'Cluster04','Cluster05','Cluster06', ...
    'Cluster07', 'Cluster08', 'Cluster09','Cluster10','Cluster11','Cluster12', ...
    'Cluster13', 'Cluster14','Cluster15','Cluster16','Cluster17','Cluster18'},savnam) ;
if fn1~= fn2
    rows = 11 ; cols = 2 ;
    figure('Position', [1 1 450*cols 40*rows]) 
    xgap = .03 ; ygap = .02 ; % percentage of total dimension
    ha = tight_subplot(rows,cols, [ygap xgap],[.01 .01],[.01 .01]) ;
    l = 0 ; r = 0; % plot counters
end

for fn = fn1:fn2
    F = []; clusters = [];
    load(fish(fn).name,'ovpmask','sumF')
    % fprintf('size of ovpmask : %d\n',size(ovpmask)) ;
    ni1 = find(ovpmask(j1,:)>= ot) ; % neurons belonging to mid
    ni2 = find(ovpmask(j2,:)>= ot) ; % neurons belonging to mid
    nni = intersect(ni1,ni2) ;
    N = length(nni) ;    % Neurons in specific mask
    cf = 0.1 ;
    fname = erase(fish(fn).name,' neurons3Dreg.mat') ;
    if isfile(cat(2,fname,' tailneuron.mat'))
        load(cat(2,fname,' tailneuron.mat'))
    elseif isfile(cat(2,dir0,'tailneuron',filesep,fname,' tailneuron.mat'))
        load(cat(2,dir0,'tailneuron',filesep,fname,' tailneuron.mat'))
    end
    [n,m] = size(F) ; % F is z-scored fluorescence of neurons (rows) as a function of time (columns)
    if max(max(F)) <= 1 
        % Z-score normalization
        zF = zscore(full(F),[],2) ;
        F = zF - min(zF,2) ; % Makes all values positive
    end    
    % if ~isfile(sprintf('%s tailangle.%s',fname,ext))            
    %     figure('Position', [1 1 600 100],'Visible','off');
    %     plot(Ts',tailneuron','Color',[0,0,0]); axis tight ; 
    %     yticks(linspace(0,90,4)) ; ylim([0 90]) ;          
    %     %title(sprintf('%s tail angle', fname)) ;       
    %     set(gca,'FontSize',10) ;
    %     set(gca,'FontName','Arial') ;        
    %     %xlabel('Time (seconds)') ; 
    %     %ylabel('deg°') ;
    %     export_fig(sprintf('%s tailangle.%s',fname,ext),'-transparent','-painters','-nocrop') ;
    %     close
    % end
    idx = clusters{tailneuron_cluster}' ; 
    idx(end) = [] ; % remove tailangle from activity matrix
    
    nrn_idx = intersect(nni,idx) ;         
    ni = length(nrn_idx) ; % tailneuron count
    % If there are tail neurons with negative tail correlation, then remove
    % cr = nan(ni,1) ;
    % for i = 1:length(nrn_idx)
    %     cr(i) = round(corr(tailneuron',F(nrn_idx(i),:)','Type','Spearman'),1) ;        
    % end
    % ii = find(cr<=0) ;
    % nrn_idx(ii) = [] ;
    % ni = length(nrn_idx) ; % revised tailneuron count    
    writematrix(fname,savnam,'Range',sprintf('A%d',fn+1)) ;
    writematrix(N,savnam,'Range',sprintf('B%d',fn+1)) ;
    writematrix(ni,savnam,'Range',sprintf('C%d',fn+1)) ;
    writematrix(round(ni/N,3),savnam,'Range',sprintf('D%d',fn+1)) ;
    
    cp = F(nrn_idx,:) ;
    scp = sum(cp(:,100:end),2) ;
    if ni > 1
        % Measure average Spearman correlation
        [I,J] = sort(scp) ;         
        cf = round(corr(tailneuron',mean(cp)','Type','Spearman'),3) ;
    elseif ni == 1
        cf = round(corr(tailneuron',cp','Type','Spearman'),3) ;
        J = 1 ;
    end
    % If there are over 500 neurons in the cluster, then split 
    % If there are non-tail neurons with tail correlation > 0.9*cp, 
    % then combine with tailneurons
    nidt = 1:n ; % All neuron indices
    nidt(idx) = [] ; % remove tailneuron indices
    nnrn_idx = intersect(nni,nidt) ;
    % for i = 1:length(nnrn_idx)
    %     cr1 = round(corr(tailneuron',F(nnrn_idx(i),:)','Type','Spearman'),3) ;
    %     if cr1 > 0.9*cf
    %         %fprintf('Neuron %d in %s added to tailneuron cluster\n',nnrn_idx(i),fname)
    %         nrn_idx = [nrn_idx nnrn_idx(i)] ;
    %     end
    % end
    % % Recalculate tail-neuron correlations and sum
    % cp = F(nrn_idx,:) ;
    % scp = sum(cp(:,100:end),2) ;
    % if ni > 1
    %     % Measure Spearman correlation
    %     [I,J] = sort(scp) ;         
    %     cf = round(corr(tailneuron',mean(cp)','Type','Spearman'),3) ;
    % elseif ni == 1
    %     cf = round(corr(tailneuron',cp','Type','Spearman'),3) ;
    %     J = 1 ;
    % end
    tna = round(sum(sumF(nrn_idx)),3) ; % tailneuron activity sum
    writematrix(tna,savnam,'Range',sprintf('E%d',fn+1)) ;            
    if ni>0
        writematrix(cf,savnam,'Range',sprintf('F%d',fn+1)) ;   
    end
        
    nn = length(nnrn_idx) ;
    if nn > 0        
        ncp = F(nnrn_idx,:) ;
        nscp = sum(ncp(:,100:end),2) ;    
        if nn > 1
            [I1,J1] = sort(nscp) ; 
            % Measure Spearman correlation
            cnf = round(corr(tailneuron',mean(ncp)','Type','Spearman'),3) ;
        elseif nn == 1
            cnf = round(corr(tailneuron',ncp','Type','Spearman'),3) ;
            J1 = 1 ;
        end
        ntna = round(sum(sumF(nnrn_idx)),3) ; % tailneuron activity sum
    else
        ntna = 0 ;
        cnf = nan ;
    end
    writematrix(ntna,savnam,'Range',sprintf('G%d',fn+1)) ;            
    writematrix(cnf,savnam,'Range',sprintf('H%d',fn+1)) ;    

    if N > 1 && prnt == 1
        % Let's print all tailneuron clusters    
        if ni > 0
            fprintf('%d of %d neurons in %s are tail-correlated\n',ni,N,fname) 
             if ni > 1
                activity = sum(cp) ;
                activity = activity./N ;
            elseif ni == 1
                activity = cp./N ;
            end
            if fn1 == fn2
                figure('Position', [1 1 600 100],'Visible','off') ;
            else  
                if contains(fname,'ctrl')
                    l = l + 1 ; % control on left
                    ik = (l-1)*2 + 1 ;
                    axes(ha(ik)) ;
                else
                    r = r + 1 ; % defective on right
                    ik = r*2 ;
                    axes(ha(ik)) ;
                end
            end
            plot(Ts',activity,'Color',[1,0,0.2]); axis tight ;
        end
        if nn > 0
            if nn > 1
                activity = sum(ncp) ;
                % only for 70 and 102, no multiplier
                activity = activity./N ;
            elseif nn == 1
                activity = ncp./N ;
            end 
            hold on ;
            if nn < 50
                plot(Ts',activity,'Color',[0.2,0,1]); axis tight ; 
            end
        end
        % ylim([0 pmx])
        % axis off
        if ni > 0 || nn > 0
            title(sprintf('%s | n=%d,%d dfF=%.4g,%.4g cc=%.4g,%.4g', fname,ni,nn,tna,ntna,cf,cnf)) ;
        end
        % if nn > 3000
        %     yticks(linspace(0,3000,4)) ; ylim([0 3000]) ;
        % elseif nn > 1000
        %     yticks(linspace(0,1500,4)) ; ylim([0 1500]) ;
        % elseif nn > 200
        %     yticks(linspace(0,450,4)) ; ylim([0 450]) ;            
        % elseif nn > 50
        %     yticks(linspace(0,150,4)) ; ylim([0 150]) ;
        % else
        %     yticks(linspace(0,40,3)) ; ylim([0,40]) ;
        % end
        %title(sprintf('%s \nTotal activity= %d of %d non-locomotor neurons',mnam,round(sum(sum(cp))),nn)) ;     
               
                   
        % yyaxis right
        % plot(Ts',tailneuron','Color',[0,0,0]); axis tight ; 
        % yticks(linspace(0,90,4)) ; ylim([0 90]) ;    
        set(gca,'FontSize',8) ;
        set(gca,'FontName','Arial') ; 
        %xlabel('Time (seconds)') ;
        %ylabel(sprintf('%s',fname)) ;
        if fn1 == fn2
            export_fig(sprintf('%s%s tail-neuron cluster Sum in %s.%s',outdir,mnam,fname,ext),'-transparent','-painters','-nocrop') ;
            close 
        end
    end
    % Let's see if any non-tail clusters are interesting
    % clid = find(nclust > 100);
    % clid = clid(~ismember(clid, tailneuron_cluster)) ;
    
    nt = zeros(1,length(nclust)) ; ci = [0.5,0,0.1] ;
    if nn >= 50
        for k = 1:length(nclust)
            nrn_idx = intersect(nni,clusters{k}) ;
            nt(k) = length(nrn_idx) ;
            if k ~= tailneuron_cluster && nt(k) > 5
                ncp = F(nrn_idx,:) ;
                activity = sum(ncp) ;                
                activity = 5*activity./N ;
    
                ci(3) = min(1,ci(3) + 0.1) ;
                if ci(3) == 1
                    ci(2) = min(1,ci(2) + 0.1) ;
                end
                plot(Ts',activity,'Color',ci); 
                % legend(sprintf('Cluster%1.2d',k))  
            end
        end
    end
    hold off ; axis tight off ; 
    if ni > 500
        ylim([0 0.2*pmx]) 
    else
        ylim([0 pmx]) 
    end
    writematrix(length(nclust),savnam,'Range',sprintf('I%d',fn+1)) ;
    writematrix(nt,savnam,'Range',sprintf('J%d',fn+1)) ;               
end
if fn1 ~= fn2
    export_fig(sprintf('%s%s cluster Sum %s.%s', ...
        outdir,mnam,datestr(now,'yymmdd_HHMM'),ext),'-transparent','-painters','-nocrop') ;
    close
end
toc
end