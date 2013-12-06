% This is the superscript that has the paths for all of the patient data.

% Paths.
cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/
setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

setenv ( 'PATHPT' , '/workdir/Patient0006/007/opt' );
pathpt = getenv ( 'PATHPT' );

%load index.txt
index = 1;

[metric] = fast_temperature_obj_fxn ( path22, pathpt, index );
%[metric,diff_Iso] = test_obj_fxn ( path22, pathpt, index );

% index = index + 1;
% csvwrite ('index.txt' , index);