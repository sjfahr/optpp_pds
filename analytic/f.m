% y = f(x) for Rosenbrock
function y = f(~)
%f1 = 10*(x(2) - x(1)^2);
%f2 = x(1) - 1;

cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/
setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

setenv ( 'PATHPT' , '/workdir/Patient0006/007/opt' );
pathpt = getenv ( 'PATHPT' );

load index.txt

[metric] =  obj_fxn_unified_optics ( path22, pathpt, index );

index = index + 1;
csvwrite ('index.txt' , index);

y =  metric;
% disp('perfusion')
% x(1)
% disp('k_0')
% x(2)
% disp('mu_a')
% x(3)
