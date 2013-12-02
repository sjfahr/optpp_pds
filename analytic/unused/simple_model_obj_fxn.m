% This is the updated Bioheat_script that should be used with DF's DAKOTA
% run.

function [metric] = simple_model_obj_fxn222 ( path22, iteration );
cd( path22);
patient_path = '/workdir/Patient0002/000/opt';
patient_opt_path = strcat( path22, patient_path);
cd( patient_opt_path);
input_param = 'optpp_pds.in.';
index = num2str(iteration);
input_filename = strcat( input_param, index);
input_filename = strcat( input_filename, '.mat' );

load(input_filename);
% aaa=strtrim(aaa);
% aaa=regexp(aaa,'\s+','split'); % Separate the label from the data into two columns.

% Write every string as a number
probe_u = str2num(probe_init);
g_anisotropy = str2num(anfact_healthy);
mu_a = str2num(mu_a_healthy);
mu_s = str2num(mu_s_healthy);
k_cond = str2num(k_0_healthy);
w_perf = str2num(w_0_healthy);
% x_disp = str2num(aaa{8}{1});
% y_disp = str2num(aaa{9}{1});
% z_disp = str2num(aaa{10}{1});
% x_rot = str2num(aaa{11}{1});
% y_rot = str2num(aaa{12}{1});
% z_rot = str2num(aaa{13}{1});

robin_co=0; %dummy var

% Load the recorded power
power_log = load ('time_then_power.csv');

% Define the domain and scaling
mod_point.x=171;  % x position on image
mod_point.y=171;
mod_point.z=1;

% Import the VTK header info
path_append = strrep ( patient_opt_path, '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/', '');
path_append = strrep ( path_append, 'opt', '');
FOV_path = strcat( '/FUS4/data2/BioTex/BrainNonMDA/processed/' , path_append );
FOV_path = strcat( FOV_path , 'laser' );
cd ( FOV_path );
FOV_import = csvimport ( 'FOV.csv' );
fov = FOV_import ( 2, : );
fov = cell2mat (fov);

matrix.x = fov(1);
matrix.y = fov(2);
matrix.z = fov(3);

spacing.x = fov(4);
spacing.y = fov(5);
spacing.z = fov(6);

FOV.x = matrix.x * spacing.x;
FOV.y = matrix.y * spacing.y;
FOV.z = matrix.z * spacing.z;

% For now, x and y scaling must be equal; z = 1
scaling.x = 1;
scaling.y = 1;
scaling.z = 1;

% Build the domain
[dom,~,~]=modeled_domain(FOV,matrix,scaling,mod_point);

% Define the source; source.n must be greater than 1 to make sense. Odd is
% better than even
source.n=5;
source.length=0.01;  %~0.033 is when n=5 is visible
source.laser=linspace((-source.length/2),(source.length/2),source.n);

% Run the Bioheat model with the unique powers
[tmap_unique]=Bioheat1D( power_log,dom,source,w_perf,k_cond,g_anisotropy,mu_a,mu_s,probe_u,robin_co);
tmap_unique=tmap_unique+37;
tmap_unique(:,:,:,1)=37;

% Make the full temperature history with the unique tmaps
% [tmap]=Build_tmap_history(tmap_unique,delta_P);

% Set delta time 'dt' for dose and run it
% dt = 2.5;
% [w,~]=ArrDose_for_1D(tmap,dt);
% w ( w > 4 ) = 4;  %Put a cap on the dose

% Change directory to load Arrhenius dose and MATLAB registration data.
dose_path = strcat( '/FUS4/data2/BioTex/BrainNonMDA/processed/' , path_append );
dose_path = strcat( dose_path , 'matlab' );
cd ( dose_path );
load 'arrheniusDose.mat'
MRTI_dose_size=size(arrheniusDose.mean);
load ( 'VOI.mat' );

% Resize the tmap_unique model into the same spacing as the MRTI
aa = imresize (tmap_unique , 1/scaling.x);
bb = aa(round(matrix.x/2),round(matrix.y/2),1,:);
[~,dd] = max (bb);

aa_size = size ( aa(:,:,1,dd) );
size_diff=[(MRTI_dose_size(1)-aa_size(1)) (MRTI_dose_size(2)-aa_size(2))];
upper_left_mod = zeros((size(aa,1)+size_diff(1)),(size(aa,2)+size_diff(2)));
upper_left_mod(1:size(aa,1),1:size(aa,2)) = aa(:,:,1,dd);

% Define the intervals that will be compared
x_range   = [ (VOI.center_in_pix.y - floor(aa_size(1)/2)) (VOI.center_in_pix.y + floor(aa_size(1)/2))];
y_range   = [ (VOI.center_in_pix.x - floor(aa_size(2)/2)) (VOI.center_in_pix.x + floor(aa_size(2)/2))];

roi_x    = [ (VOI.center_in_pix.y - 60) (VOI.center_in_pix.y + 60) ]; %For model 57 deg C isotherms.
roi_y    = [ (VOI.center_in_pix.x - 60) (VOI.center_in_pix.x + 60) ];

roi_x_MRTI   = [ (VOI.center_in_pix.y - 20) (VOI.center_in_pix.y + 20) ];  %For MRTI dose.
roi_y_MRTI   = [ (VOI.center_in_pix.x - 20) (VOI.center_in_pix.x + 20) ];

x_range = round ( x_range ); % Round everything after calculating them.
y_range = round ( y_range );
roi_x = round ( roi_x );
roi_y = round ( roi_y );
roi_x_MRTI = round ( roi_x_MRTI );
roi_y_MRTI = round ( roi_y_MRTI );

% Register the model data to the MRTI
matched_mod = zeros (MRTI_dose_size(1), MRTI_dose_size(2));
matched_mod ( x_range(1):x_range(2), y_range(1):y_range(2) ,:,:) = upper_left_mod( 1:aa_size(1) , 1:aa_size(2) );

model_Iso = (matched_mod >57 ) ; % 57 deg C isotherm

model_Iso ( 1:roi_x(1), : ) = 0;
model_Iso ( roi_x(2):end,: )  = 0;
model_Iso ( :, 1:roi_y(1) ) = 0;
model_Iso ( :,roi_y(2):end ) = 0;

MRTI_Iso = (arrheniusDose.mean >1 ) ;

MRTI_Iso ( 1:roi_x_MRTI(1), : ) = 0;
MRTI_Iso ( roi_x_MRTI(2):end, : )  = 0;
MRTI_Iso ( :, 1:roi_y_MRTI(1) ) = 0;
MRTI_Iso ( :,roi_y_MRTI(2):end ) = 0;

diff_Iso= model_Iso - MRTI_Iso;

metric = abs(sum(sum(diff_Iso)));

cd (path22);
% output_param = 'optpp_pds.out.';
% index = num2str(iteration);
% output_filename = strcat( output_param, index);
% 
% csvwrite ( output_filename, metric);

end