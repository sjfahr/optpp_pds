% y = f(x) for Rosenbrock
function y = f(~)
%f1 = 10*(x(2) - x(1)^2);
%f2 = x(1) - 1;

cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/
setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

index = load ( 'index.txt' );
patient_index = load ( 'patient_index.txt' );

Patient_Paths = importdata( 'patient_paths.txt' );

setenv ( 'PATHPT' , char( Patient_Paths( patient_index ) ) ); % char is to make the cell string become an array string
pathpt = getenv ( 'PATHPT' );

[metric] =  fast_temperature_obj_fxn ( path22, pathpt, index );

index = index + 1;
csvwrite ('index.txt' , index);

y =  metric;
% disp('perfusion')
% x(1)
% disp('k_0')
% x(2)
% disp('mu_a')
% x(3)
