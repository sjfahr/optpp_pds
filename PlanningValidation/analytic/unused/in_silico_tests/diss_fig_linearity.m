[hot_true] = Human_GPU_choice_sym ( 15,spacing,scaling,mod_point,source,w_perf,k_cond,g_anisotropy,mu_eff_list,probe_u,robin_co,c_blood,choice);

hot_linearity = all_opt_fig .* 15 - (15 - 1).*no_pwr_fig;

diff = (hot_true-hot_linearity);
max(diff(:))
norm(diff,2)

figure; imagesc(no_pwr_fig);
set(findobj('type','axes'),'fontsize',14);
xlabel('Pixel number in ROI (Unity)')
ylabel('Pixel number in ROI (Unity)')
h = colorbar;
set(findobj('type','axes'),'fontsize',14);
ylabel(h,strcat('Temperature (',sprintf('%cC', char(176)),')'));

figure; imagesc(all_opt_fig);
set(findobj('type','axes'),'fontsize',14);
xlabel('Pixel number in ROI (Unity)')
ylabel('Pixel number in ROI (Unity)')
h = colorbar;
set(findobj('type','axes'),'fontsize',14);
ylabel(h,strcat('Temperature (',sprintf('%cC', char(176)),')'));

figure; imagesc(hot_linearity);
set(findobj('type','axes'),'fontsize',14);
xlabel('Pixel number in ROI (Unity)')
ylabel('Pixel number in ROI (Unity)')
h = colorbar;
set(findobj('type','axes'),'fontsize',14);
ylabel(h,strcat('Temperature (',sprintf('%cC', char(176)),')'));

figure; imagesc(diff);
set(findobj('type','axes'),'fontsize',14);
xlabel('Pixel number in ROI (Unity)')
ylabel('Pixel number in ROI (Unity)')
h = colorbar;
set(findobj('type','axes'),'fontsize',14);
ylabel(h,strcat('Temperature (',sprintf('%cC', char(176)),')'));