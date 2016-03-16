clear
close all
clc
tic
setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

data_filename = 'datasummaryL2_10sourceNewton50.txt';  % Name the datasummary file

opttype = 'bestfit50' ;
datasummary = dlmread(data_filename,',',1,0);
datasummary(any(isnan(datasummary), 2), 7) = 1;
num_studies = size(datasummary,1);

opt.paths = cell (num_studies,2);
    
Study_paths{1,1} = strcat( 'Study00',num2str(datasummary(1,1)));
Study_paths{1,2} = strcat( '0',num2str(datasummary(1,2)));


% From mu_eff_data, find the matching study's(ies') mu_eff value(s)
total = cell(2,7);
mat_string = ['workdir/',Study_paths{1,1},'/',Study_paths{1,2},'/opt/optpp_pds.',opttype,'.in.1.mat'];
load (mat_string);

n_sources = [1e+0 1e+1 1.5e+1 2.5e+1 5e+1 7.5e+1 1e+2 2e+2 1e+3 1e+4 1e+5];
%n_sources = [1e+0 1.2e+1 1e+1 1e+2];
total{1,1}{1} = [Study_paths{1,1},'/',Study_paths{1,2}];
total{2,1}{1} = [Study_paths{1,1},'/',Study_paths{1,2}];
total{1,1}{2} = 200;
total{2,1}{2} = 2000;

xx= inputdatavars.voi(2)-inputdatavars.voi(1)+1;
yy= inputdatavars.voi(4)-inputdatavars.voi(3)+1;
l_n_sources = length(n_sources);
inputdatavars.powerhistory = '[[36, 42, 49, 142, 160],[0.0, 3.0, 0.0, 15 0.0]]';
for ii = 1:2
    mu = total{ii,1}{2};
    temp_fields = zeros(yy,xx,l_n_sources);
    for jj = 1:l_n_sources
%         % Display run information
%         disp('Start ')
%         disp(strcat (num2str(ii),' of ', num2str(num_studies)))
%         fprintf('iter %s \n', total{ii,1});
        
        % Do global optimization
        [temp_fields(:,:,jj)] = temperature_obj_fxn_GPU_tmap ( inputdatavars, n_sources(jj), mu );
        
    end
    total{ii,2}=temp_fields;
end
clear ii jj
cd ../../../MATLAB/Tests/in_silico/
save ('in_silico_test22.mat','total');
cd (path22);

for ii = 1:2
    total{ii,3} = zeros(l_n_sources,3);
    total{ii,4} = zeros(l_n_sources,15);
    total{ii,5} = total{ii,4};
    total{ii,6} = total{ii,4};
    total{ii,7} = total{ii,4};

    for jj = 1:l_n_sources
        
        tmap_iter = total{ii,2}(:,:,jj);
        tmap_gold = total{ii,2}(:,:,end);
        t_diff = tmap_gold - tmap_iter;

        total{ii,3}(jj,1)= n_sources(jj);
        total{ii,3}(jj,2)= ( norm ( t_diff , 2 ) )^2; % L2 norm
        total{ii,3}(jj,3)= max(max( total{ii,2}(:,:,jj)  ) ) - max( tmap_gold(:));  % 
        
        for kk = 1:15
            model_deg_threshold = tmap_iter >= ( 50 + kk);
            [mod_row, mod_column] = find( model_deg_threshold ==1);
            mod_list = [(mod_row .* inputdatavars.spacing(1)) (mod_column .* inputdatavars.spacing(2))];
            gold_deg_threshold = tmap_gold >= ( 50 + kk);
            [gold_row, gold_column] = find( gold_deg_threshold ==1);
            gold_list = [(gold_row .* inputdatavars.spacing(1)) (gold_column .* inputdatavars.spacing(2))];
            n_model = sum( sum( model_deg_threshold ));
            n_gold =  sum( sum( gold_deg_threshold  ));
            intersection = model_deg_threshold + gold_deg_threshold;
            intersection = intersection > 1;
            n_intersection = sum( sum( intersection ));
            total{ii,4}(jj,kk) = n_model - n_intersection; % False positive count
            total{ii,5}(jj,kk) = n_gold - n_intersection; % False negative count
            total{ii,6}(jj,kk) = 2 * n_intersection / ( n_model + n_gold );  % DSC
            total{ii,7}(jj,kk) = HausdorffDist( mod_list, gold_list);
        end
        
    end
end

cd ../../../MATLAB/Tests/in_silico/
save ('in_silico_test15W.mat','total');
cd (path22);

clear ii jj
for ii = 1:l_n_sources
    figure;imagesc(total{1,2}(:,:,ii),[37,max(max(total{1,2}(:,:,1)))])
    set(findobj('type','axes'),'fontsize',14);
    xlabel('Pixel number in ROI (Unity)')
    ylabel('Pixel number in ROI (Unity)')
    h = colorbar;
    set(findobj('type','axes'),'fontsize',14);
    ylabel(h,strcat('Temperature (',sprintf('%cC', char(176)),')'));
end
clear ii h
for ii = 1:l_n_sources-1
    figure;imagesc(total{1,2}(:,:,10)-total{1,2}(:,:,ii))
    set(findobj('type','axes'),'fontsize',14);
    xlabel('Pixel number in ROI (Unity)')
    ylabel('Pixel number in ROI (Unity)')
    h = colorbar;
    set(findobj('type','axes'),'fontsize',14);
    ylabel(h,strcat('Temperature difference (',sprintf('%cC', char(176)),')'));
end
clear ii h
aa=total{1,3}(1:end-1,:);
fig=figure('units','normalized','position',[.1 .3 .18 .45]);
[AX,H1,H2] = plotyy(log10(aa(:,1)),log10(aa(:,2)),log10(aa(:,1)),log10(abs(aa(:,3))));
set(findobj('type','axes'),'fontsize',15);
set(H1,'linestyle','--','marker','x'); set(H2,'linestyle','-.','marker','o');
set([H1;H2],'MarkerSize',15);
xlabel('Log_{10} (M Sources) (Unity)')
ylabel(AX(1),'Log_{10} (L_2 norm) (Unity)')
ylabel(AX(2),strcat('Log_{10} [abs(Max Diff)] [Log_{10} (',sprintf('%cC', char(176)),')]'))
print(fig,'mu200c','-dpng');

aa=total{2,3}(1:end-1,:);
fig=figure('units','normalized','position',[.1 .3 .18 .45]);
%fig=figure('units','normalized','position',[.1 .3 .22 .55]);
[AX,H1,H2] = plotyy(log10(aa(:,1)),log10(aa(:,2)),log10(aa(:,1)),log10(abs(aa(:,3))));
set(findobj('type','axes'),'fontsize',15);
set(H1,'linestyle','--','marker','x'); set(H2,'linestyle','-.','marker','o');
set([H1;H2],'MarkerSize',15);
xlabel('Log_{10} (M Sources) (Unity)')
ylabel(AX(1),'Log_{10} (L_2 norm) (Unity)')
ylabel(AX(2), strcat('Log_{10} [abs(Max Diff)] [Log_{10} (',sprintf('%cC', char(176)),')]'))
print(fig,'mu2000c','-dpng');

for ii =1:2
    for jj = 1:10
        
        tmap_iter = total{ii,2}(:,:,jj);
        tmap_gold = total{ii,2}(:,:,end);
        t_diff = tmap_gold - tmap_iter;
        
        total{ii,3}(jj,1)= n_sources(jj);
        total{ii,3}(jj,2)= ( norm ( t_diff , 2 ) )^2; % L2 norm
        total{ii,3}(jj,3)= max(max( total{ii,2}(:,:,jj)  ) ) - max( tmap_gold(:));  %
        total{ii,3}(jj,4)= ( norm ( t_diff , inf) )^2; %L_inf norm
        
        
    end
end