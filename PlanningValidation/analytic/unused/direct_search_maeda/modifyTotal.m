% This script adds two sections to 'total' by presenting the optimal
% modified objective function
clear ii jj
close all
%clear

setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

cd ../../../MATLAB/Tests/direct_search/libraries

choice=1;
toss_choice = 0;
metric_choice = 1;  % This is only relevant if choice == 5. 1 is DSC, 2 is L2
index_range_choice = 1;
index_range = [1 401];

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
    
    %load ('GPU_dict_perf_mu_global_400.mat');
    load ('GPU_dict_perf_mu_global_400');
    
elseif choice ==5  % random
     
    load ('GPU_dict_perf_mu_rand.mat');
    
    
end

cd(path22);

%total(1,:)=[];

if toss_choice == 1
    ix       =find(~cellfun(@isempty,regexp(total(:,1),'0497'))==1); % Absolutely should be excluded  % good to exclude (median)
    
    %ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0378'))==1); %Strongly suggest exclusion %  good to keep (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0476'))==1); % Strongly suggest exclusion  % good to exclude (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0435'))==1); % very probably suggest exclusion: susceptibility artifact; new  % good to exclude (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0440'))==1); % very probably suggest exclusion: susceptibility artifact;new  % good to exclude (median)
    
    
    
    %ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0436'))==1); %good to keep (ParaView)  % good to keep (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0466'))==1); % Very probably suggest exclusion % good to exclude (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0468'))==1); % Very probably suggest exclusion % good to exclude (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0471'))==1); % Strongly suggest exclusion % good to exclude (median)
    
    %ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0417'))==1); %Very probably suggest exclusion % good to keep (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0409'))==1); % Absolutely should be excluded % good to exclude (median)
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0415'))==1); % Absolutely should be excluded % good to exclude (median)
    total(ix,:) = [];
    
    % ix=find(~cellfun(@isempty,regexp(total(:,1),'0457'))==1);
    % total(ix,:) = [];
    
    % ix=find(~cellfun(@isempty,regexp(total(:,1),'0457'))==1);
    % total(ix,:) = [];
end

pt_length = size(total,1);
n_length = size (total{2,2},1);
total = [total cell(pt_length,2)];
total{1,12} = 'L2norm with Dice penalty';
total{1,13} = 'Optimal L2norm with Dice penalty';

for ii=2:pt_length
    total{ii,12} = ( total{ii,2}(:,3) + 1 ./ (total{ii,3}(:,7) + 1E-7) ) ./ 2; % 1/(dice + 1E-7)
end
clear ii

if index_range_choice == 1

    for ii=2:pt_length
        total{ii,13} = zeros(1,3);
        [total{ii,13}(1) , index] = min (total{ii,12}(index_range(1):index_range(2))); % record optimal L2+Dice penalty
        total{ii,13}(2) = total{ii,2}(index,1); % record value that produces optimal L2+Dice penalty
        total{ii,13}(3) = index;
    end
    
else
    for ii=2:pt_length
        total{ii,13} = zeros(1,3);
        [total{ii,13}(1) , index] = min (total{ii,12}(:)); % record optimal L2+Dice penalty
        total{ii,13}(2) = total{ii,2}(index,1); % record value that produces optimal L2+Dice penalty
        total{ii,13}(3) = index;
    end
    
end

cd ../../../MATLAB/Tests/direct_search/libraries
save ('GPU_dict_mu_mod400.mat','total')
cd(path22);
