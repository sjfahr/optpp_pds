% This is the superscript that has the paths for all of the patient data.

% Paths. Be in the folder you want to process. E.g.:
% '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0002/000/opt'
setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

load index.txt

[metric] = simple_model_obj_fxn ( path22, index );

index = index + 1;
csvwrite ('index.txt' , index);