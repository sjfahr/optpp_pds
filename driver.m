close all
clear all

%load data
load simulatedData.mat

npixel = 8
roi = [2,7;
       2,7;]
InitialGuess = GenerateInitialGuess (MRData,npixel,roi);

LowerBound = zeros(size(InitialGuess)); % pixel loc, pixel loc, params: amplitude, shape, scale, delay
UpperBound = zeros(size(InitialGuess)); % pixel loc, pixel loc, params: amplitude, shape, scale, delay

LowerBound(:,:,1) = 0; % amplitide: 0-10
LowerBound(:,:,2) = 2 ; % shape : 2-7
LowerBound(:,:,3) = 1 ; % scale : 1-3
LowerBound(:,:,4) = 0; % delay: 0-40
UpperBound(:,:,1) = 10; % amplitide: 0-10
UpperBound(:,:,2) = 7 ; % shape : 2-7
UpperBound(:,:,3) = 3 ; % scale : 1-3
UpperBound(:,:,4) = 40; % delay: 0-20

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
options = optimset('jacobian','off','MaxIter',10000000000,'MaxFunEvals',10000000000);
[x,resnorm,residual,exitflag,output,lambda] = lsqnonlin( @(x) dynProjKernel(x,MRData,npixel,roi),InitialGuess(:),LowerBound(:),UpperBound(:),options);

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

%% Generate data based on guess
%
%theta = 0;
%for proj = 1:N,
%    theta = theta + 111.246;
%    im = squeeze(Bp(:,:,proj));
%    Rp(:,proj) = radon(im,theta);
%end
%figure;imshow(Rp,[])

%%% Compare errors in the sinogram
%error = sqrt(sum(sum((Rp - R).*conj(Rp - R))))
%rmsError = dynProjKernel(Ap(:,:,1),Ap(:,:,2),Ap(:,:,3),Ap(:,:,4),R)
%
%figure;imshow(abs([Rp ; R]),[]);

%% modify guess to more closely match "true" sinogram


