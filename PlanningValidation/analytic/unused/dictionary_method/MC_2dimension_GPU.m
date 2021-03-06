% This script finds the best mu_eff for the different studies.
function MC_2dimension_GPU (choice); % Monte Carlo_2dimension_GPU
%choice = 1; % 1 = mu; 2 = perf; 3 = cond;
tic

cd /mnt/FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
data_filename = 'datasummaryL2_10sourceNewton50.txt';  % Name the datasummary file

opttype = 'bestfit50';
datasummary = dlmread(data_filename,',',1,0);
datasummary(any(isnan(datasummary), 2), 7) = 1;
num_studies = size(datasummary,1);
for ii = 1:num_studies
    
    Study_paths{ii,1} = strcat( 'Study00',num2str(datasummary(ii,1)));
    Study_paths{ii,2} = strcat( '0',num2str(datasummary(ii,2)));
    
end
clear ii

[~, ~, max_phys_sz]=display_inputvars(1:30);

% choice = 1 is mu; 2 is perf; 3 is cond;
input_path = cell(1,2);
input_path{1,1} = Study_paths{1,1};
input_path{1,2} = Study_paths{1,2};
[ all_opt_fig, no_pwr_fig,sim_dim, summary ] = Generate_MC_library ( max_phys_sz, choice, input_path, opttype );

cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/direct_search/libraries

if choice == 4
    
    save ('MC_mu.mat','all_opt_fig','no_pwr_fig','sim_dim','summary','-v7.3');
    
% elseif choice == 2
%     
%     save ('all_opt_perf.mat','all_opt_fig','no_pwr_fig','sim_dim','summary','-v7.3');
%     
% elseif choice == 3
%     
%     save ('all_opt_cond.mat','all_opt_fig','no_pwr_fig','sim_dim','summary','-v7.3');
%     
% elseif choice ==4
%     
%     save ('all_opt_perf_mu.mat', 'all_opt_fig','no_pwr_fig','sim_dim','summary','-v7.3');
%     
% elseif choice ==5
%     
%     save ('rand_opt_perf_mu.mat', 'all_opt_fig','no_pwr_fig','sim_dim','summary','-v7.3');
end
cd /mnt/FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
end