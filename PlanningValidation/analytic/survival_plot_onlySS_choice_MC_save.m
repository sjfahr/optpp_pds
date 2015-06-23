%
function [] = survival_plot_onlySS_choice_MC_save (dice_values, naive_values, dice_opt,naive_tag);

%close all

sizes.var = size(dice_values,1);
sizes.dice = size(dice_values,2);

thresholds = linspace ( 0.0, 1, 10001);
passes = cell(sizes.var,sizes.dice);
naive_pass = passes;
L2_mean_pass =passes;
L2_median_pass = passes;
dice_mean_pass = passes;
dice_median_pass = passes;
for ii = 1:sizes.var
    
    for jj = 1:sizes.dice
        passes{ii,jj} = zeros(10001,1);
        naive_pass{ii,jj} = passes{ii,jj};
%         L2_mean_pass{ii,jj} = passes{ii,jj};
%         L2_median_pass{ii,jj} = passes{ii,jj};
%         dice_mean_pass{ii,jj} = passes{ii,jj};
%         dice_median_pass{ii,jj} = passes{ii,jj};
        

        
    end
end
clear ii jj

passes_opt = zeros(10001,1);
pass_naive_opt = passes_opt;
passes_LOOCV = passes_opt;
% all_L2_mean_pass = passes_opt;
% all_L2_median_pass = passes_opt;
% all_dice_mean_pass = passes_opt;
% all_dice_median_pass = passes_opt;
for kk = 1:10001
    
    passes_opt (kk) = sum( dice_opt > thresholds(kk) );
    pass_naive_opt (kk) = sum ( naive_values > thresholds(kk) );
    passes_LOOCV (kk) = sum(dice_values > thresholds(kk) );
%     all_L2_mean_pass (kk) = sum( best.all.L2_mean.val > thresholds(kk) );
%     all_L2_median_pass (kk) = sum( best.all.L2_median.val > thresholds(kk) );
%     all_dice_mean_pass (kk) = sum( best.all.dice_mean.val > thresholds(kk) );
%     all_dice_median_pass (kk) = sum( best.all.dice_median.val > thresholds(kk) );
end
clear kk

% figure; h_title = title( 'DSC performance for optimization versus literature and single best guesses'); hold all;
% [h1] = plot (thresholds, [passes_opt pass_naive_opt all_L2_mean_pass all_L2_median_pass all_dice_mean_pass all_dice_median_pass]);
% legend( h1, 'Optimization', strcat(['Literature '], num2str(naive_tag(1))), 'L_2 mean','L_2 median','dice mean','dice median') ;
% legend('-DynamicLegend', 'Location','southwest');hold off;
cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/display_performance/survival_plots
fig=figure('units','normalized','position',[.1 .3 .22 .52]);
hold all;
[h1] = plot (thresholds, [passes_opt pass_naive_opt],'LineWidth',5);
set(findobj('type','axes'),'fontsize',15);
xlim( [0.4 1]);
xlabel('DSC (Unity)');
ylabel('Number of passing datasets (Unity)');
set(findobj('type','axes'),'fontsize',15);
print(fig,'opt_survival20','-dpng');


%         xlabel('\mu_{eff}   [ m^{-1} ]'); ylabel('\rho [ kg/(m^3 s) ]'); set(findobj('type','axes'),'fontsize',15);
%         print(fig,'opt_map','-dpng');
%         var_opt(ii,:)
%         Study_paths{ii,2}

fig=figure('units','normalized','position',[.1 .3 .22 .52]);
hold all;
[h1] = plot (thresholds, [passes_LOOCV pass_naive_opt],'LineWidth',5);
set(findobj('type','axes'),'fontsize',15);
xlim( [0.4 1]);
xlabel('DSC (Unity)');
ylabel('Number of passing datasets (Unity)');
set(findobj('type','axes'),'fontsize',15);
print(fig,'LOOCV_survival20','-dpng');

%legend( h1, 'Optimization', ['Literature ', num2str(naive_tag(1)),' m^{-1}' ]) ;
%legend('-DynamicLegend', 'Location','southwest');hold off;

% kk = 2;
% 
% for ii = 1:sizes.var
%     
% 
% 
%         
%         figure;
%         
%         h_title = title( ['DSC thresholds of ']);
%         hold on
%         for jj = 1:sizes.dice
%             
%          
%                 hold all
%                 origtitle = get(h_title,'String');
%                 set(h_title, 'String', strcat(origtitle, [', DSC > 0']));
%                 set(findobj('type','axes'),'fontsize',13);
%                 
%                
%                     
%                     if jj == 1
%                         plot(thresholds, passes{ii,jj},'LineWidth',5);
%                         legend_string = ['DSC > 0'];
%                         hold all
% 
%                         
%                         if naive_tag ==1
% 
%                             plot (thresholds, naive_pass{ii,jj}, 'DisplayName',legend_string, 'LineWidth',5);
%                         end
%                         
%                     else
% 
%                         set(findobj('type','axes'),'fontsize',13);
%                         plot(thresholds, passes{ii,jj},'DisplayName',legend_string,'LineWidth',5);
% 
%                     end
%                 
%                 
%                 
%             
%         end
%         
%         hold off
%         kk=kk+1;
%         
%     
% end

end
