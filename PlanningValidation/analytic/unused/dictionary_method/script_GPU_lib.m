close all
clear
clc
choice = 4; %Dataset selection
quick_choice=2; %shortened processing

GPU_generate_lib(choice);

metric_generation(choice, quick_choice);

%MC_2dimension_GPU (choice);

% for ii = 2:3
%     GPU_generate_lib(ii);
%     
%     metric_generation(ii);
% end