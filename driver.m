close all
clear all

options = optimset('jacobian','off');
MRData=rand(10);
LowerBound   = [-1 -1];
UpperBound   = [10 10];
InitialGuess = [0,1];
%If FUN is parameterized, you can use anonymous functions to capture the
%problem-dependent parameters. Suppose you want to solve the non-linear
%least squares problem given in the function myfun, which is
%parameterized by its second argument c. Here myfun is a MATLAB file
%function such as
%
%        function F = myfun(x,c)
%        F = [ 2*x(1) - exp(c*x(1))
%              -x(1) - exp(c*x(2))
%              x(1) - x(2) ];
%
%To solve the least squares problem for a specific value of c, first
%assign the value to c. Then create a one-argument anonymous function
%that captures that value of c and calls myfun with two arguments.
%Finally, pass this anonymous function to LSQNONLIN:
%
%        c = -1; % define parameter first
%        x = lsqnonlin(@(x) myfun(x,c),[1;1])
[x,resnorm,residual,exitflag,output,lambda] = lsqnonlin(@(x) ForwardProject(x,MRData),InitialGuess,LowerBound,UpperBound,options)

switch exitflag
 case 1
   disp('Function converged to a solution x.')
 case 2 
   disp('Change in x was less than the specified tolerance.')
 case 3  
   disp('Change in the residual was less than the specified tolerance.')
 case 4  
   disp('Magnitude of search direction was smaller than the specified tolerance.')
 case 0  
   disp('Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.MaxFunEvals.')
 case -1 
   disp('Output function terminated the algorithm.')
 case -2 
   disp('Problem is infeasible: the bounds lb and ub are inconsistent.')
 case -4 
   disp('Line search could not sufficiently decrease the residual along the current search direction.')
 otherwise
   disp('unknown ')
end
x
%MRData(1)
