% y = f(x) for Rosenbrock
function y = f(x)
%f1 = 10*(x(2) - x(1)^2);
%f2 = x(1) - 1;

setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

index = load ( 'index.txt' );

[metric] = simple_model_obj_fxn ( path22, index );

index = index + 1;
csvwrite ('index.txt' , index);

y =  metric;
disp('perfusion')
x(1)
disp('k_0')
x(2)
disp('mu_a')
x(3)
