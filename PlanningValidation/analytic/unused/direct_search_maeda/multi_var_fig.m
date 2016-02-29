clear ii jj
close all
clear

setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

cd ../../../MATLAB/Tests/direct_search

choice=1;
if choice == 1   % mu
    
    load ('GPU_global_mu2.mat');
    
elseif choice == 2  % perf
    
    load ('GPU_global_perf2.mat');
    
elseif choice == 3   % cond
    
    load ('GPU_global_cond2.mat');
    
end

cd(path22);

total(1,:)=[];

aa = cell2mat( total(:,8));

if choice == 1   % mu
    
    index=find( (aa(:,2)>700)==1);
    
elseif choice == 2  % perf
    
    index=find( (aa(:,2)<26)==1);
    
elseif choice == 3   % cond
    
    index=find( (aa(:,2)>0.4)==1);
    
end
sz = size(aa);
if length (sz) == 2
    
    bb = zeros( length(index), sz(2));
    
    
    for jj=1:length(index)
        ii = index(jj);
        figure(ii);
        [AX h1 h2] = plotyy(total{ii,2}(:,1),total{ii,2}(:,3),total{ii,2}(:,1),total{ii,3}(:,7));
        legend([h1;h2],'L_2','DSC');
        bb(ii,:)=aa(index,:);
    end
    figure;
    bb(index,:)=[];
    hist(bb(:,2));
    var=Descriptive_statistics(bb(:,2))
    DSC=Descriptive_statistics_LOOCV(bb(:,1))
end

if length (sz) == 3
    
    bb = zeros( length(index), sz(2),sz(3));
    
    
    for jj=1:length(index)
        ii = index(jj);
        figure(ii);
        [AX h1 h2] = plotyy(total{ii,2}(:,1),total{ii,2}(:,3),total{ii,2}(:,1),total{ii,3}(:,7));
        legend([h1;h2],'L_2','DSC');
        bb(ii,:,:)=aa(index,:,:);
    end
    figure;
    aa(index,:)=[];
    hist(aa(:,2));
    var=Descriptive_statistics(aa(:,2))
    DSC=Descriptive_statistics_LOOCV(aa(:,1))
end
