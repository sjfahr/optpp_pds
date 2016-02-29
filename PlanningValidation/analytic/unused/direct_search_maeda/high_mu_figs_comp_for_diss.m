clear ii jj
close all
%clear

setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

% cd ../../../MATLAB/Tests/direct_search/
cd ../../../MATLAB/Tests/direct_search/libraries

choice=1;
toss_choice = 1;
metric_choice = 1;  % This is only relevant if choice == 5. 1 is DSC, 2 is L2
index_range_choice = 1;
index_range = [1 401];

if choice == 1   % mu
    
    %load ('GPU_global_mu2.mat');
    if index_range_choice ==1
        load('GPU_dict_mu_mod400.mat');
    else
        load ('GPU_dict_mu.mat');
    end
    
elseif choice == 2  % perf
    
    %load ('GPU_global_perf2.mat');
    load ('GPU_dict_perf.mat');
    
elseif choice == 3   % cond
    
    %load ('GPU_global_cond2.mat');
    load ('GPU_dict_cond.mat');
    
elseif choice ==4 
    
    %load ('GPU_dict_perf_mu_global_400.mat');
    load ('GPU_dict_perf_mu_global_400');
    
elseif choice ==5  % random
     
    load ('GPU_dict_perf_mu_rand.mat');
    
    
end

cd(path22);

total(1,:)=[];

if toss_choice == 1
    ix       =find(~cellfun(@isempty,regexp(total(:,1),'0497'))==1); % Absolutely should be excluded  % good to exclude (median)
    
    %ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0378'))==1); %Strongly suggest exclusion %  good to keep (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0476'))==1); % Strongly suggest exclusion  % good to exclude (median)
    
    %ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0435'))==1); % very probably suggest exclusion: susceptibility artifact; new  % good to exclude (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0440'))==1); % very probably suggest exclusion: susceptibility artifact;new  % good to exclude (median)
    
    
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0436'))==1); %good to keep (ParaView)  % good to keep (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0466'))==1); % Very probably suggest exclusion % good to exclude (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0468'))==1); % Very probably suggest exclusion % good to exclude (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0471'))==1); % Strongly suggest exclusion % good to exclude (median)
    
    %ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0417'))==1); %Very probably suggest exclusion % good to keep (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0409'))==1); % Absolutely should be excluded % good to exclude (median)
    
    %ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0415'))==1); % Absolutely should be excluded % good to exclude (median)
    total(ix,:) = [];
    
    % ix=find(~cellfun(@isempty,regexp(total(:,1),'0457'))==1);
    % total(ix,:) = [];
    
    % ix=find(~cellfun(@isempty,regexp(total(:,1),'0457'))==1);
    % total(ix,:) = [];
end

aa = cell2mat( total(:,8));

if choice == 1   % mu
    
    index=find( (aa(:,2)>-1)==1);
    
elseif choice == 2  % perf
    
    index=find( (aa(:,2)<26)==1);
    
elseif choice == 3   % cond
    
    index=find( (aa(:,2)>0.4)==1);
    
elseif choice == 4
    index=find( (aa(:,2)>-1)==1);
elseif choice ==5
    index=find( (aa(:,2)>-1)==1);
    
end

if choice ==5
    
    w_lim = [3.25 16];
    %w_lim = [3 16.5];
    w_pix = 100;
    w_array = linspace( w_lim(1),w_lim(2),w_pix);
    mu_lim = [55 2500];
    %mu_lim = [50 3000];
    mu_pix = 200;
    mu_array = linspace( mu_lim(1),mu_lim(2),mu_pix);
    %[paraXq, paraYq] = meshgrid (w_lim(1): 0.25 : w_lim(2), mu_lim(1): 50: mu_lim(2) );
    [paraXq, paraYq] = meshgrid ( w_array, mu_array);
    
%     bb=hist2(total{ii,2}(:,1),total{ii,2}(:,2),200);  
%     figure;imagesc(bb); xlabel('perf');ylabel('mu');
    
    for jj=1:length(index)
        ii =index(jj);
        
        if metric_choice ==1
            [x, y, z]=griddata(total{ii,2}(:,1),total{ii,2}(:,2),total{ii,3}(:,7),paraYq,paraXq);
            figure; contourf(x,y,z);colorbar; caxis([0 0.9]); title(['DSC ' total{ii,1}]); xlabel('mu_{eff}   [ m^{-1} ]'); ylabel('perf [ kg/(m^3 s) ]');
        elseif metric_choice==2
            [x, y, z]=griddata(total{ii,2}(:,1),total{ii,2}(:,2),total{ii,2}(:,4),paraYq,paraXq);
            figure; contourf(x,y,z);colorbar; title(['L_2 ' total{ii,1}]); xlabel('mu_{eff}   [ m^{-1} ]'); ylabel('perf [ kg/(m^3 s) ]');
        end
        
        
    end
    
elseif choice==4

    
    [paraXq, paraYq] = meshgrid ( unique(summary.w_perf), unique(summary.mu));
    
    for jj=1:length(index)
        ii =index(jj);
        dice_iter = squeeze( total{ii,3});
        
        if metric_choice ==1
            % [x, y, z]=griddata(total{ii,2}(:,1),total{ii,2}(:,2),total{ii,3}(:,7),paraYq,paraXq);
            % figure; contourf(x,y,z);colorbar; caxis([0 0.9]); title(['DSC ' total{ii,1}]); xlabel('mu_{eff}   [ m^{-1} ]'); ylabel('perf [ kg/(m^3 s) ]');
            [x, y, z]=griddata(total{ii,2}(:,1),total{ii,2}(:,2),dice_iter,paraYq,paraXq);
            figure; contourf(x,y,z);colorbar; caxis([0 0.9]); title(['DSC ' total{ii,1}]); xlabel('mu_{eff}   [ m^{-1} ]'); ylabel('perf [ kg/(m^3 s) ]');
        elseif metric_choice==2
            [x, y, z]=griddata(total{ii,2}(:,1),total{ii,2}(:,2),total{ii,2}(:,4),paraYq,paraXq);
            figure; contourf(x,y,z);colorbar; title(['L_2 ' total{ii,1}]); xlabel('mu_{eff}   [ m^{-1} ]'); ylabel('perf [ kg/(m^3 s) ]');
        end
        
        
    end
    
else   % This has been modded for comparing DAKOTA's opt with global
    
    AX=zeros(length(index),2);
    h1=zeros(length(index),1);
    h2 = h1;
    
    
    for jj=1:length(index)
        ii = index(jj);
        figure(ii);
        [AX(ii,:) h1(ii) h2(ii)] = plotyy(total{ii,2}(index_range(1):index_range(2),1),total{ii,12}(index_range(1):index_range(2)),...
            total{ii,2}(index_range(1):index_range(2),1),total{ii,3}(index_range(1):index_range(2),7));
        legend([h1(ii);h2(ii)],'Mod obj fxn','DSC');
    end
    figure;
    %aa(index,:)=[];
    hist(aa(:,2));
    var=Descriptive_statistics(aa(:,2))
    DSC=Descriptive_statistics_LOOCV(aa(:,1))
    
end

clear ii jj
cd /FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation

% List the studies and number of optimization realizations
% 0389; 0402; 0451; 0453
path = cell(4,1);
path{1} = './workdir/Study0018/0402/opt/';
path{2} = './workdir/Study0018/0389/opt/';
path{3} = './workdir/Study0026/0453/opt/';
path{4} = './workdir/Study0026/0451/opt/';

num_real = zeros(4,1);
num_real(1) = 13;
num_real(2) = 14;
num_real(3) = 13;
num_real(4) = 8;

real_mu_runs = cell(4,1);

for ii=1:4
    real_mu_runs{ii}=zeros(num_real(ii),1);
    for jj = 1:num_real(ii)
        load (strcat(path{ii},'optpp_pds.bestfit50.in.',num2str(jj),'.mat'));
        real_mu_runs{ii}(jj) = str2num(inputdatavars.cv.mu_eff_healthy);
        
    end
end

for ii = [7 8 16 18] % 0402; 0389; 0453; 0451

        fig=figure(ii);
        figure(ii);
        set(AX(ii,1),'YLim',[100 500]);
        set(AX(ii,1),'YTick',[100:100:500]);
        set(AX(ii,2),'YLim',[0 1]);
        set(AX(ii,2),'YTick',[0:.25:1]);
        %print(fig,'compDAKOTA_global','-dpng');
        hold on
        plot( 180, total{ii,12}(181), 'bo')
        
        if ii == 7
            
            cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/display_performance/0402
        elseif ii == 8
            cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/display_performance/0378
        elseif ii ==16
            cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/display_performance/0453
        elseif ii == 18
            cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/display_performance/0451
        end
        print(fig,'compDAKOTA_global','-dpng');
end
cd (path22);
