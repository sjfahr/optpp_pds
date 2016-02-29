clear

cd /FUS4/data2/sjfahrenholtz/MATLAB/Tests

load brute_sense_Aug_14

figure;plot(total_path{1,2}(:,1),total_path{1,3}(:,7),total_path{1,2}(:,1),total_path{1,3}(:,1),total_path{1,2}(:,1),total_path{1,3}(:,15))

cd /FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation

[max_dice,max_index]=max(total_path{1,3}(:,7));

Study_paths=cell(1,2);
Study_paths{1}='Study0026';
Study_paths{2}='0455';

aaa=load('/FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation/workdir/Study0026/0455/opt/optpp_pds.bestfit4.in.2.mat');
aaa.inputdatavars.cv.mu_eff_healthy = num2str(total_path{1,2}(max_index,1));
[L2norm, dice, tmap_model, MRTI_crop] =  temperature_obj_fxn ( aaa.inputdatavars, 10 );

model_deg = tmap_model > 57;
MRTI_deg = 2*(MRTI_crop > 57);
n_model = sum(sum( model_deg ));
n_MRTI = sum(sum( MRTI_deg ));
intersection = model_deg + MRTI_deg;

figure; imagesc(tmap_model,[30 80]);
figure; imagesc(MRTI_crop, [30 80]);
figure; imagesc(intersection)


%[L2norm, dice, tmap_model, MRTI_crop,mat_struct,quality] = check_opt_GoodBadUgly( Study_paths, 'bestfit4', datasummary(:,3));

