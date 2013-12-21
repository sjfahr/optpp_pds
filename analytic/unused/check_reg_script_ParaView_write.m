% This is the superscript that allows you to input the VOI information.

% Paths. Be in the folder you want to process. E.g.:
% '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0002/000/opt'



cd /FUS4/data2/BioTex/BrainNonMDA/processed
setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
setenv ( 'WORKDIR' , pwd) ;
work_dir = getenv ( 'WORKDIR' );

patient_index = 1;
Patient_Paths = importdata( 'patient_paths.txt' );
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );

VOI_pre.x = [ 83 123] ;
VOI_pre.y = [ 108 148 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 103 128 0 ];
VOI_pre.time = 60;

% time = 143;
% crop_big = 25;
% crop_small = 5;

%[ center, VOI ] = check_reg ( path22 ,time,crop_big,crop_small );
%[ center, VOI ] = check_reg_find_center ( path22 ,VOI_pre );

[ center, VOI ] = check_reg_Paraview_write (  VOI_pre, path22, work_dir, pathpt );