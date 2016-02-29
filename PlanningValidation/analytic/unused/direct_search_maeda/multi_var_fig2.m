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

Lnorms= cell2mat( total(:,7));
Lnorms= Lnorms(:,2);
nn=3;
L1norms= Lnorms(1: nn : end);
L2norms= Lnorms(2: nn : end);
Linfnorms= Lnorms(3: nn : end);

DSC = cell2mat( total(:,8));
HD = cell2mat( total(:,9));

MI = cell2mat( total(:,10));
MI = MI(:,2);
nn = 2;
MI_temp = MI(1: nn : end);
MI_isotherm = MI(2: nn : end);

Fpixel = cell2mat( total(:,11));

for ii=1
    figure;
    [h1] = plot(total{ii,2}(:,1),total{ii,2}(:,2:4) );
    title('Distance norms');
    legend(h1,'L_1','L_2','L_{inf}');
    
    figure;
    [h1] = plot(total{ii,2}(:,1), [total{ii,3}(:,1) total{ii,3}(:,7) total{ii,3}(:,15)]);
    title('DSC');
    legend(h1,'51 ^{o}C','57 ^{o}C','65 ^{o}C');
    
    figure;
    [h1] = plot(total{ii,2}(:,1), [total{ii,4}(:,1) total{ii,4}(:,7) total{ii,4}(:,15)]);
    ylim([0 0.03]);
    title('Hausdorf Distance');
    legend(h1,'51 ^{o}C','57 ^{o}C','65 ^{o}C');
    
    figure;
    [AX h1 h2] = plotyy(total{ii,2}(:,1), [total{ii,5}(:,1) total{ii,5}(:,7) total{ii,5}(:,15)],total{ii,2}(:,1),total{ii,2}(:,5));
    title('Mutual Information');
    legend([h1;h2],'51 ^{o}C','57 ^{o}C','65 ^{o}C','Temperature MI');
    
    figure;
    [h1] = plot(total{ii,2}(:,1), [total{ii,6}(:,1,3) total{ii,6}(:,7,3) total{ii,6}(:,15,3)]);
    title('Total false pixel count (total false = FN + FP)');
    legend([h1;h2],'51 ^{o}C','57 ^{o}C','65 ^{o}C');
    
    figure;
    [AX h1 h2] = plotyy(total{ii,2}(:,1),total{ii,2}(:,3),total{ii,2}(:,1),total{ii,3}(:,7));
    title('My preferred comparison: L_2 and DSC');
    legend([h1;h2],'L_2','DSC');
    %bb(ii,:)=aa(index,:);
    
end


figure;
hist(L1norms);
title('Optimal mu_{eff} value for L_1 norm');
L1_var=Descriptive_statistics(L1norms)

figure;
hist(L2norms);
title('Optimal mu_{eff} value for L_2 norm');
L2_var=Descriptive_statistics(L2norms)

figure;
hist(Linfnorms);
title('Optimal mu_{eff} value for L_{inf} norm');
Linf_var=Descriptive_statistics(Linfnorms)

%bb(index,:)=[];
figure;
hist(DSC(:,2));
title('Optimal mu_{eff} value for DSC, 57 ^{o}C');
DSC_var=Descriptive_statistics(DSC(:,2))
DSC_stat=Descriptive_statistics_LOOCV(DSC(:,1))

figure;
hist(HD(:,2));
title('Optimal mu_{eff} value for HD, 57 ^{o}C');
HD_var=Descriptive_statistics(HD(:,2))

figure;
hist(MI_isotherm);
title('Optimal mu_{eff} value for MI, 57 ^{o}C');
MI_isotherm_var=Descriptive_statistics(MI_isotherm)

figure;
hist(MI_temp);
title('Optimal mu_{eff} value for Temperature field');
MI_temp_var=Descriptive_statistics(MI_temp)

figure;
hist(Fpixel(:,2));
title('Optimal mu_{eff} value for false pixels (FN+FP) for , 57 ^{o}C');
Fpixel_var=Descriptive_statistics(Fpixel(:,2))
