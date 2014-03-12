%  x = vector of input parameters
%  y = vector of output projection data
function y = ForwardProject(x,data)
f1 = 10*(x(2) - x(1)^2);
f2 = x(1) - 1;
y = [f1,f2];
%disp('mrdata test')
%data(1)
