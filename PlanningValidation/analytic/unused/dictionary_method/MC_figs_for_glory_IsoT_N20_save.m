clear ii jj
close all
clear

setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

% cd ../../../MATLAB/Tests/direct_search/
cd ../../../MATLAB/Tests/direct_search/libraries

choice=4;
toss_choice = 1;
metric_choice = 1;  % This is only relevant if choice == 5. 1 is DSC, 2 is L2
if choice == 1   % mu
    
    %load ('GPU_global_mu2.mat');
    load ('GPU_dict_mu.mat');
    
elseif choice == 2  % perf
    
    %load ('GPU_global_perf2.mat');
    load ('GPU_dict_perf.mat');
    
elseif choice == 3   % cond
    
    %load ('GPU_global_cond2.mat');
    load ('GPU_dict_cond.mat');
    
elseif choice ==4
    
    %load ('GPU_dict_perf_mu_global_400'); 
    %load ('GPU_dict_perf_mu_global_400_all_metricRedo');
    %load ('GPU_dict_perf_mu_global_400_all_metric');
    %load ('GPU_dict_perf_mu_global_400_all_metric_3');
    %load ('opt_perf_mu_400_long_diss');
    load ('long_diss_ArrDose22');
elseif choice ==5  % random
    
    load ('GPU_dict_perf_mu_rand.mat');
    
    
end

cd(path22);

total(1,:)=[];

if toss_choice == 1
    ix       =find(~cellfun(@isempty,regexp(total(:,1),'0497'))==1); % Absolutely should be excluded  % good to exclude (median)
    
    %ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0378'))==1); %Strongly suggest exclusion %  good to keep (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0476'))==1); %Not positive;; % Strongly suggest exclusion  % good to exclude (median)
    
%     ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0385'))==1); %temporary hack to check if registration is what needs to change
%     
%     ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0495'))==1); %temporary hack to check if registration is what needs to change
    
    % for STM 2015 con   %ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0435'))==1); % Not positive % very probably suggest exclusion: susceptibility artifact; new  % good to exclude (median);
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0435'))==1); % Not positive % very probably suggest exclusion: susceptibility artifact; new  % good to exclude (median);
    
    %for STM 2015 con   %ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0440'))==1); % very probably suggest exclusion: susceptibility artifact;new  % good to exclude (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0440'))==1); % very probably suggest exclusion: susceptibility artifact;new  % good to exclude (median)
    
    %ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0436'))==1); %good to keep (ParaView)  % good to keep (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0466'))==1); % Very probably suggest exclusion % good to exclude (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0468'))==1); % Very probably suggest exclusion % good to exclude (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0471'))==1); % Strongly suggest exclusion % good to exclude (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0447'))==1); % Not positive% Very probably suggest exclusion % good to exclude (median)  
    
    %ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0417'))==1); %Very probably suggest exclusion % good to keep (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0409'))==1); % Absolutely should be excluded % good to exclude (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0415'))==1); % Not positive% Absolutely should be excluded % good to exclude (median) 
    total(ix,:) = [];
    
    % ix=find(~cellfun(@isempty,regexp(total(:,1),'0457'))==1);
    % total(ix,:) = [];
end

aa = cell2mat( total(:,8)); %This sets it seek the best 

if choice == 1   % mu
    
    index=find( (aa(:,2)>-1)==1);
    
elseif choice == 2  % perf
    
    index=find( (aa(:,2)<26)==1);
    
elseif choice == 3   % cond
    
    index=find( (aa(:,2)>0.4)==1);
    
    
elseif choice ==5 || choice ==4
    index=find( (aa(:,2)>-1)==1);
    
end

if choice ==5 ||choice==4
    
    if choice==4
        w_array = unique( summary.w_perf);
        mu_array = unique( summary.mu);
    else
        w_lim = [3.25 16];
        %w_lim = [3 16.5];
        w_pix = 250;
        w_array = linspace( w_lim(1),w_lim(2),w_pix);
        mu_lim = [55 400];
        %mu_lim = [50 3000];
        mu_pix = 250;
        mu_array = linspace( mu_lim(1),mu_lim(2),mu_pix);
    end
    
    
    %[paraXq, paraYq] = meshgrid (w_lim(1): 0.25 : w_lim(2), mu_lim(1): 50: mu_lim(2) );
    [paraXq, paraYq] = meshgrid ( w_array, mu_array);
    
    %     bb=hist2(total{ii,2}(:,1),total{ii,2}(:,2),200);
    %     figure;imagesc(bb); xlabel('perf');ylabel('mu');
    grid_sz = size (paraXq);
    Xx = zeros(grid_sz(1), grid_sz(2), length(index));
    Yy = Xx;
    obj_fxn = Xx;
    for jj=1:length(index)
        ii =index(jj);
        %dice = squeeze(total{ii,3});
        dice = total{ii,3}(:,7);
        %dice = total{ii,12}(:,7);
        if metric_choice ==1
            %[Xx(:,:,jj), Yy(:,:,jj), obj_fxn(:,:,jj)]=griddata(total{ii,2}(:,1),total{ii,2}(:,2),total{ii,3}(:,7),paraYq,paraXq);
            [Xx(:,:,jj), Yy(:,:,jj), obj_fxn(:,:,jj)]=griddata(total{ii,2}(:,1),total{ii,2}(:,2),dice,paraYq,paraXq);
            %figure; contourf(x,y,z);colorbar; caxis([0 0.9]); title(['DSC ' total{ii,1}]); xlabel('mu_{eff}   [ m^{-1} ]'); ylabel('perf [ kg/(m^3 s) ]');
            %         elseif metric_choice==2
            %             [x, y, z]=griddata(total{ii,2}(:,1),total{ii,2}(:,2),total{ii,2}(:,4),paraYq,paraXq);
            %             figure; contourf(x,y,z);colorbar; title(['L_2 ' total{ii,1}]); xlabel('mu_{eff}   [ m^{-1} ]'); ylabel('perf [ kg/(m^3 s) ]');
        end
        
        
        
    end
    clear ii jj
    
    cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/display_performance/survival_plots/N20_IsoT
    opt_mean = mean(obj_fxn,3);
    fig=figure('units','normalized','position',[.1 .3 .264 .624]);
    contourf(Xx(:,:,1),Yy(:,:,1),opt_mean);caxis([0 0.9]); colorbar; %title(['Global Mean DSC']);
    h = colorbar;
    set(findobj('type','axes'),'fontsize',15);
    ylabel(h,strcat('DSC (Unity)'));
    %set(findobj('type','axes'),'fontsize',15);
    xlabel('\it\mu_{eff} \rm(m^{-1})'); ylabel('\it\omega \rm( kg/(m^3 s) )');
    print(fig,'mean_global','-dpng');
    %naive_point = [
    
    opt_median = median(obj_fxn,3);
    fig=figure('units','normalized','position',[.1 .3 .264 .624]);
    contourf(Xx(:,:,1),Yy(:,:,1),opt_median);caxis([0 0.9]); colorbar; %title(['Global Median DSC']);
    xlabel('mu_{eff}   [ m^{-1} ]'); ylabel('perf [ kg/(m^3 s) ]'); set(findobj('type','axes'),'fontsize',15);
    h = colorbar;
    set(findobj('type','axes'),'fontsize',15);
    ylabel(h,strcat('DSC (Unity)'));
    %set(findobj('type','axes'),'fontsize',15);
    xlabel('\it\mu_{eff} \rm(m^{-1})'); ylabel('\it\omega \rm( kg/(m^3 s) )');
    print(fig,'median_global','-dpng');
    hold on
    naive_spot=[180 6];
    scatter( naive_spot(1),naive_spot(2) ,400,'d','markerfacecolor','b','markeredgecolor','g',...
        'LineWidth',2);
    hold off
    print(fig,'median_global_naive20','-dpng');
    
    opt_std = std(obj_fxn,0,3);
    fig=figure('units','normalized','position',[.1 .3 .264 .624]);
    contourf(Xx(:,:,1),Yy(:,:,1),opt_std); caxis([0 0.2]); colorbar; %title(['Global St. Dev. of DSC']);
    xlabel('mu_{eff}   [ m^{-1} ]'); ylabel('perf [ kg/(m^3 s) ]'); set(findobj('type','axes'),'fontsize',15);
    h = colorbar;
    set(findobj('type','axes'),'fontsize',15);
    ylabel(h,strcat('DSC (Unity)'));
    %set(findobj('type','axes'),'fontsize',15);
    xlabel('\it\mu_{eff} \rm(m^{-1})'); ylabel('\it\omega \rm( kg/(m^3 s) )');
    print(fig,'std_global','-dpng');
    hold on
    naive_spot=[180 6];
    scatter( naive_spot(1),naive_spot(2) ,400,'d','markerfacecolor','b','markeredgecolor','g',...
        'LineWidth',2);
    hold off
    print(fig,'std_global_naive20','-dpng');
    
    opt_pass=obj_fxn;        
    opt_pass(opt_pass< 0.7)=0;
    opt_pass(opt_pass>=0.7)=1;
    opt_pass_sum = sum(opt_pass,3);
    max_pass = max(opt_pass_sum(:));
    fig=figure('units','normalized','position',[.1 .3 .264 .624]);
    contourf(Xx(:,:,1),Yy(:,:,1),opt_pass_sum,[0 2 4 6 8 10 12 14 16 max_pass] );%colorbar; title(['Map of passing datasets']);
    h = colorbar;
    set(findobj('type','axes'),'fontsize',15);
    ylabel(h,strcat('Datasets > 0.7 (Unity)'));
    %set(findobj('type','axes'),'fontsize',15);
    xlabel('\it\mu_{eff} \rm(m^{-1})'); ylabel('\it\omega \rm( kg/(m^3 s) )');
    print(fig,'pass_global','-dpng');
    hold on
    naive_spot=[180 6];
    scatter( naive_spot(1),naive_spot(2) ,400,'d','markerfacecolor','b','markeredgecolor','g',...
        'LineWidth',2);
    hold off
    print(fig,'pass_global_naive20','-dpng');
    cd /mnt/FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
    
    clear ii jj
% %     jj=1;
%     opt_union=[1 1 1];
%     intersect=ones(size(obj_fxn,1),size(obj_fxn,2));
% %     while max( max( sum(opt_union,3))) < length(index)& jj<=80
%         thresh = 0.7;% - 0.01 *jj;
%         aa=cell2mat(total(:,8));
%         aa=aa(:,1);
%         aa=thresh*ones(size(aa));
%         opt_union = obj_fxn;
%         for ii = 1:length(index)
%             iter = opt_union(:,:,ii);
%             iter( iter< aa(ii))=0;
%             iter( iter>=aa(ii))=1;
%             opt_union(:,:,ii) = iter;
%             tmp(:,:)=obj_fxn(:,:,ii);
%             tmp=heaviside(tmp-thresh).*tmp;
%             tmp(tmp>0)=1;
%             intersect=intersect.*tmp;
%             
%             %         opt_union(opt_union<aa(ii))=0;
%             %         opt_union(opt_union>aa(ii))=1;
%             %         opt_pass22 = sum(opt_union,3);
%         end
% %         jj=jj+1;
% %     end
%         opt_pass22 = sum(opt_union,3);
%     figure; contourf(Xx(:,:,1),Yy(:,:,1),opt_pass22);colorbar; title(['Map of passing datasets']);
%     xlabel('mu_{eff}   [ m^{-1} ]'); ylabel('perf [ kg/(m^3 s) ]'); set(findobj('type','axes'),'fontsize',15);
        
    obj_fxn_sz = size(obj_fxn);
    LOOCV_mean = zeros( obj_fxn_sz(1),obj_fxn_sz(2),obj_fxn_sz(3) );
    LOOCV_median = LOOCV_mean;
    LOOCV_std = LOOCV_mean;
    mean_ix1 = zeros( length(index),1);
    mean_ix2 = mean_ix1;
    mean_max = mean_ix1;
    median_ix1 = mean_ix1;
    median_ix2 = mean_ix1;
    median_max = mean_ix1;
    pass_ix1 = mean_ix1;
    pass_ix2 = mean_ix1;
    pass_max = mean_ix1;
    II_mean = mean_ix1;
    II_median = mean_ix1;
    II_pass = mean_ix1;
    no_global_mu = mean_ix1;
    no_global_perf = mean_ix1;
    no_global_ix1 = mean_ix1;
    no_global_ix2 = mean_ix1;
    prc_max = zeros (length(index),4);
    prc_ix1=prc_max;
    prc_ix2=prc_max;
    II_prc=prc_max;
    prc_top_max = zeros (length(index),4);
    prc_top_ix1=prc_max;
    prc_top_ix2=prc_max;
    II_top_prc=prc_max;

    w_perf_list = unique(summary.w_perf);
    mu_list = unique(summary.mu);
    for jj=1:length(index)
        ii = index(jj);
        
        obj_fxn_iter = obj_fxn;
        obj_fxn_iter( :,:,ii ) = [];
        obj_fxn_iter_pass = obj_fxn_iter;
        obj_fxn_iter_pass( obj_fxn_iter_pass < 0.7)=0;
        obj_fxn_iter_pass( obj_fxn_iter_pass >=0.7)=1;
        obj_fxn_iter_pass= sum(obj_fxn_iter_pass,3);
        rep_max = max(obj_fxn_iter_pass(:));
        
        no_global_iter = aa;
        no_global_iter(ii,:)=[];
        no_global_mu_iter = no_global_iter(:,2)/2;
        no_global_mu(ii) = 2 * round(mean(no_global_mu_iter));
        no_global_perf_iter = no_global_iter(:,3)*4;
        no_global_perf(ii) = round(mean(no_global_perf_iter)) / 4;
        
        
%         figure; contourf(Xx(:,:,1),Yy(:,:,1),obj_fxn_iter_pass );colorbar; title(['Map of passing datasets']);
%         xlabel('mu_{eff}   [ m^{-1} ]'); ylabel('perf [ kg/(m^3 s) ]'); set(findobj('type','axes'),'fontsize',15);
        
        obj_fxn_iter_pass( obj_fxn_iter_pass < rep_max-1) = 0;
        obj_fxn_iter_pass( obj_fxn_iter_pass >=rep_max-1) = 1;
        rep_iter_pass = repmat( obj_fxn_iter_pass, [1 1 size(obj_fxn_iter,3)]);
        obj_fxn_top_pass = obj_fxn_iter .* rep_iter_pass;
        
        obj_fxn_top_prc = prctile(obj_fxn_top_pass,[20 30 40 50 60 70 80 90],3);
        for kk=1:size(obj_fxn_top_prc,3)
            prc_top_iter = obj_fxn_top_prc(:,:,kk);
            [prc_top_max(ii,kk), II_top_prc(ii,kk)] = max( prc_top_iter(:));
            [prc_top_ix1(ii,kk), prc_top_ix2(ii,kk)]=ind2sub( size(prc_top_iter),II_top_prc(ii,kk));
        end
        
        
        LOOCV_mean = mean(obj_fxn_iter,3);
        LOOCV_median = median( obj_fxn_iter,3);
        
        obj_fxn_iter_pass = obj_fxn_iter;
        obj_fxn_iter_pass( obj_fxn_iter_pass < 0.7 ) = 0;
        obj_fxn_iter_pass( obj_fxn_iter_pass >=0.7 ) = 1;
        LOOCV_pass = sum( obj_fxn_iter_pass, 3);
        
        obj_fxn_prc = prctile(obj_fxn_iter,[20 30 40 50 60 70 80 90],3);
        for kk=1:size(obj_fxn_prc,3)
            prc_iter = obj_fxn_prc(:,:,kk);
            [prc_max(ii,kk), II_prc(ii,kk)] = max( prc_iter(:));
            [prc_ix1(ii,kk), prc_ix2(ii,kk)]=ind2sub( size(prc_iter),II_prc(ii,kk));
        end
        
        [mean_max(ii), II_mean(ii)] = max( LOOCV_mean(:));
        [mean_ix1(ii), mean_ix2(ii)] = ind2sub( size(LOOCV_mean), II_mean(ii) );
        
        [median_max(ii), II_median(ii)] = max( LOOCV_median(:));
        [median_ix1(ii), median_ix2(ii)] = ind2sub( size(LOOCV_median), II_median(ii) );
        
        [pass_max(ii), II_pass(ii)] = max( LOOCV_pass(:));
        [pass_ix1(ii), pass_ix2(ii)] = ind2sub( size(LOOCV_pass), II_pass(ii) );
        
        no_global_ix1(ii) = find( mu_list == no_global_mu(ii));
        no_global_ix2(ii) = find( w_perf_list == no_global_perf(ii));
        
        
    end
    LOOCV_mean_pre = [ mean_max mean_ix1 mean_ix2 II_mean];      % max value; mu ix; w ix; column ix;
    LOOCV_median_pre = [ median_max median_ix1 median_ix2 II_median];
    LOOCV_pass_pre = [ pass_max pass_ix1 pass_ix2 II_pass];
    
    
%     LOOCV_mean_post = zeros (length(index),1);
%     LOOCV_median_post = LOOCV_mean_post;
%     LOOCV_pass_post = LOOCV_mean_post;
    LOOCV_mean_post = zeros (length(index),1);
    LOOCV_median_post = LOOCV_mean_post;
    LOOCV_pass_post = LOOCV_mean_post;
    naive_pass = LOOCV_mean_post;
    best_eyeball_norm = LOOCV_mean_post;
    eyeball_norm20=LOOCV_mean_post;
    no_global_LOOCV = LOOCV_mean_post;
    mean_spot = LOOCV_mean_post;
    no_global_spot=LOOCV_mean_post;
    prc_runs = zeros( length(index),4);
    prc_LOOCV = prc_runs;
    prc_top_runs = prc_runs;
    prc_top_LOOCV =prc_runs;
    
    naive_ix = find( total{ii,2}(:,1) == 180 & total{ii,2}(:,2)==6);
    eyeball_ix =  find( total{ii,2}(:,1) == 212 & total{ii,2}(:,2)==10.25);
    eyeball_ix20 =  find( total{ii,2}(:,1) == 250 & total{ii,2}(:,2)==12);
    for jj = 1:length(index)
        ii=index(jj);        
        mean_spot(ii) = find( total{ii,2}(:,1) == mu_array(LOOCV_mean_pre(ii,2)) & total{ii,2}(:,2)==w_array(LOOCV_mean_pre(ii,3)));
        median_spot = find( total{ii,2}(:,1) ==mu_array(LOOCV_median_pre(ii,2)) & total{ii,2}(:,2)==w_array(LOOCV_median_pre(ii,3)));
        pass_spot = find(total{ii,2}(:,1)==mu_array(LOOCV_pass_pre(ii,2)) & total{ii,2}(:,2)==w_array(LOOCV_pass_pre(ii,3)));
        no_global_spot(ii) = find( total{ii,2}(:,1) == mu_array(no_global_ix1(ii)) & total{ii,2}(:,2)== w_array(no_global_ix2(ii)) );
        
        %         obj_fxn_iter( :,:,ii ) = [];
        %         LOOCV_mean = mean(obj_fxn_iter,3);
        %         LOOCV_median = median( obj_fxn_iter,3);
        %
        %         obj_fxn_iter_pass = obj_fxn_iter;
        %         obj_fxn_iter_pass( obj_fxn_iter_pass < 0.7 ) = 0;
        %         obj_fxn_iter_pass( obj_fxn_iter_pass >=0.7 ) = 1;
        %         LOOCV_pass = sum( obj_fxn_iter_pass, 3);
        
        %         LOOCV_mean_post(ii) = max( max( LOOCV_mean));
        %         LOOCV_median_post(ii) = max( max(LOOCV_median));
        %         LOOCV_pass_post(ii) = max(max(LOOCV_pass));
        
%         LOOCV_mean_post(ii) = total{ii,3}(II_mean(ii));
%         LOOCV_median_post(ii) = total{ii,3}(II_median(ii));
%         LOOCV_pass_post(ii) = total{ii,3}(II_pass(ii));
        
%         LOOCV_mean_post(ii) = total{ii,3}(mean_spot);
%         LOOCV_median_post(ii) = total{ii,3}(median_spot);
%         LOOCV_pass_post(ii) = total{ii,3}(pass_spot);
%         naive_pass(ii) = total{ii,3}(naive_ix);
%         best_eyeball_norm(ii)= total{ii,3}(eyeball_ix);
%         eyeball_norm22(ii) = total{ii,3}(eyeball_ix22);
        
%         LOOCV_mean_post(ii) = total{ii,12}(mean_spot(ii),7);
%         LOOCV_median_post(ii) = total{ii,12}(median_spot,7);
%         LOOCV_pass_post(ii) = total{ii,12}(pass_spot,7);
%         naive_pass(ii) = total{ii,12}(naive_ix,7);
%         best_eyeball_norm(ii)= total{ii,12}(eyeball_ix,7);
%         eyeball_norm20(ii) = total{ii,12}(eyeball_ix20,7);
%         no_global_LOOCV(ii) = total{ii,12}(no_global_spot(ii),7);

        LOOCV_mean_post(ii) = total{ii,3}(mean_spot(ii),7);
        LOOCV_median_post(ii) = total{ii,3}(median_spot,7);
        LOOCV_pass_post(ii) = total{ii,3}(pass_spot,7);
        naive_pass(ii) = total{ii,3}(naive_ix,7);
        best_eyeball_norm(ii)= total{ii,3}(eyeball_ix,7);
        eyeball_norm20(ii) = total{ii,3}(eyeball_ix20,7);
        no_global_LOOCV(ii) = total{ii,3}(no_global_spot(ii),7);
        
        for kk = 1:size(obj_fxn_prc,3)
%             prc_runs(ii,kk) = find( total{ii,2}(:,1) == mu_array(prc_ix1(ii,kk)) & total{ii,2}(:,2)==w_array(prc_ix2(ii,kk)));
%             prc_LOOCV(ii,kk) = total{ii,3}(prc_runs(ii,kk));
%             prc_top_runs(ii,kk) = find( total{ii,2}(:,1) == mu_array(prc_top_ix1(ii,kk)) & total{ii,2}(:,2)==w_array(prc_top_ix2(ii,kk)));
%             prc_top_LOOCV(ii,kk) = total{ii,3}(prc_top_runs(ii,kk));
            
            prc_runs(ii,kk) = find( total{ii,2}(:,1) == mu_array(prc_ix1(ii,kk)) & total{ii,2}(:,2)==w_array(prc_ix2(ii,kk)));
            prc_LOOCV(ii,kk) = total{ii,2}(prc_runs(ii,kk),7);
            prc_top_runs(ii,kk) = find( total{ii,2}(:,1) == mu_array(prc_top_ix1(ii,kk)) & total{ii,2}(:,2)==w_array(prc_top_ix2(ii,kk)));
            prc_top_LOOCV(ii,kk) = total{ii,2}(prc_top_runs(ii,kk),7);
            
        end
        
        
    end
%     LOOCV_mean_post_stats = Descriptive_statistics_LOOCV(LOOCV_mean_post);
%     LOOCV_median_post_stats= Descriptive_statistics_LOOCV(LOOCV_median_post);
%     LOOCV_pass_post_stats = total{ii,3}(II_pass(ii));
    LOOCV_mean_post_stats = Descriptive_statistics_LOOCV(LOOCV_mean_post);
    LOOCV_median_post_stats= Descriptive_statistics_LOOCV(LOOCV_median_post);
    LOOCV_pass_post_stats = Descriptive_statistics_LOOCV(LOOCV_pass_post);
    naive_stats = Descriptive_statistics_LOOCV(naive_pass);
    best_eyeball_stats = Descriptive_statistics_LOOCV(best_eyeball_norm);
    eyeball_stat20=Descriptive_statistics_LOOCV(eyeball_norm20);
    no_global_LOOCV_stats = Descriptive_statistics_LOOCV(no_global_LOOCV);
    
    survival_plot_onlySS_choice_MC_IsoT_save (LOOCV_mean_post, naive_pass, aa(:,1), 1,20);
%    survival_plot_compareAll ( naive, SSgrad, FEM, SSglobal )
%     prc_LOOCV60 = Descriptive_statistics_LOOCV(prc_LOOCV(:,1));
%     prc_LOOCV70 = Descriptive_statistics_LOOCV(prc_LOOCV(:,2));
%     prc_LOOCV80 = Descriptive_statistics_LOOCV(prc_LOOCV(:,3));
%     prc_LOOCV90 = Descriptive_statistics_LOOCV(prc_LOOCV(:,4));
else
    fig=figure('units','normalized','position',[.1 .3 .264 .624]);
    contourf(Xx(:,:,1),Yy(:,:,1),opt_mean);caxis([0 0.9]); colorbar; %title(['Global Mean DSC']);
    h = colorbar;
    set(findobj('type','axes'),'fontsize',15);
    ylabel(h,strcat('DSC (Unity)'));
    %set(findobj('type','axes'),'fontsize',15);
    xlabel('\it\mu_{eff} \rm(m^{-1})'); ylabel('\it\omega \rm( kg/(m^3 s) )');
    
    
    
    
    
    
    for jj=1:length(index)
        ii = index(jj);
        figure(ii);
        [AX h1 h2] = plotyy(total{ii,2}(:,1),total{ii,2}(:,3),total{ii,2}(:,1),total{ii,2}(:,7));
        legend([h1;h2],'L_2','DSC');
    end
    figure;
    %aa(index,:)=[];
    hist(aa(:,2));
    var=Descriptive_statistics(aa(:,2))
    DSC=Descriptive_statistics_LOOCV(aa(:,1))
    
end


cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/MATLAB_file_exchange
% invprctile(aa(:,1),[0.5 0.6 0.7 0.8 0.85 0.9 0.95 0.975])
invprctile(LOOCV_mean_post,[0.5 0.6 0.7 0.8 0.85 0.9 0.95 0.975])
invprctile(naive_pass,[0.5 0.6 0.7 0.8 0.85 0.9 0.95 0.975])
cd /mnt/FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation

%Optimal points
cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/display_performance/survival_plots/N20_IsoT
aafit=aa(:,2:3);
aafit([1 3],:)=[]; % Throw out the outliers
fig=figure('units','normalized','position',[.1 .3 .264 .624]);
contourf(Xx(:,:,1),Yy(:,:,1),opt_mean);caxis([0 0.9]); colorbar; %title(['Global Mean DSC']);
h = colorbar;
set(findobj('type','axes'),'fontsize',15);
ylabel(h,strcat('DSC (Unity)'));
xlabel('\it\mu_{eff} \rm(m^{-1})'); ylabel('\it\omega \rm( kg/(m^3 s) )');
hold on
scatter(aafit(:,1),aafit(:,2),150,'MarkerFaceColor','b','MarkerEdgeColor','g',...
    'LineWidth',2);
hh=lsline;
ylim([0.5 16.5])
xlim([50 400])
set(hh(1),'color','b','LineWidth',3.5);
%set(hh(2),'LineWidth',0.0001);
%scatter(aa([1 3],2),aa([1 3],3),150,'markerfacecolor','g','markeredgecolor','g');
scatter(aa([1 3],2),aa([1 3],3),150,'markerfacecolor','g','markeredgecolor','g'); % label the outliers as green
hold off
print(fig,'global_scatter20','-dpng');

%Unique points run during LOOCV
aafit = aa(:,2:3);
fig=figure('units','normalized','position',[.1 .3 .264 .624]);
contourf(Xx(:,:,1),Yy(:,:,1),opt_mean);caxis([0 0.9]); colorbar; %title(['Global Mean DSC']);
h = colorbar;
set(findobj('type','axes'),'fontsize',15);
ylabel(h,strcat('DSC (Unity)'));
xlabel('\it\mu_{eff} \rm(m^{-1})'); ylabel('\it\omega \rm( kg/(m^3 s) )');
hold on
unique_mean_spot = unique(mean_spot);
scatter(total{1,2}(unique_mean_spot,1),total{1,2}(unique_mean_spot,2),150,'markerfacecolor','b','markeredgecolor','g',...
    'LineWidth',2);
bbfit = [ round(mean(aafit(:,1))) round( mean(aafit(:,2)).*4)./4];
indiv_opt_ix = find( total{ii,2}(:,1) == bbfit(1) & total{ii,2}(:,2)==bbfit(2));
hold off
print(fig,'global_LOOCV_choice20','-dpng');

% Mean map with naive point
fig=figure('units','normalized','position',[.1 .3 .264 .624]);
contourf(Xx(:,:,1),Yy(:,:,1),opt_mean);caxis([0 0.9]); colorbar; %title(['Global Mean DSC']);
h = colorbar;
set(findobj('type','axes'),'fontsize',15);
ylabel(h,strcat('DSC (Unity)'));
xlabel('\it\mu_{eff} \rm(m^{-1})'); ylabel('\it\omega \rm( kg/(m^3 s) )');
hold on
naive_spot=[180 6];
scatter( naive_spot(1),naive_spot(2) ,400,'d','markerfacecolor','b','markeredgecolor','g',...
    'LineWidth',2);
hold off
print(fig,'mean_global_naive20','-dpng');

cd /mnt/FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation