close all
clear all
clc

%% initialize
N = 100; % number of projections
%x = linspace(-30,120,N); % not sure if need negative
x = 1:N;


%% setup theoretical phantom
A = zeros(8,8,4); % pixel loc, pixel loc, params: amplitude, shape1, shape2, delay
B = zeros(8,8,length(x));   % full times series
for I = 2:7,
    for J = 2:7,
        A(I,J,1) = 10*rand; % amplitide: 0-10
        A(I,J,2) = 2 + 5*rand; % shape 1: 2-7
        A(I,J,3) = 1 + 2*rand; % shape 2: 1-3
        A(I,J,4) = 40*rand; % delay: 0-20
        B(I,J,:) = A(I,J,1)*gampdf(x-A(I,J,4),A(I,J,2),A(I,J,3)); % bring them into time series
    end
end
A(5:6,6:7,1) = 20*A(5:6,6:7,1); % change model amplitide to represent tumor
B(5:6,6:7,:) = 20*B(5:6,6:7,:); % represent tumor

%% Generate synthetic data
theta = 0;
for proj = 1:N,
    theta = theta + 111.246; % increment based on golden ratio
    im = squeeze(B(:,:,proj));
    R(:,proj) = radon(im,theta); % calculate synthetic acquired sinogram
end
figure;imshow(R,[])

%% recon the images using regular backprojection (optional)

% proj = 1:1:N;
% psi = (proj-1)*111.246;
% imr = iradon(R,psi); % produce the image
% %figure;imshow(imr,[])

%% make a guess to starting gamma
gammaGuess = sum(R,1);
%figure;plot(x,gammaGuess);
%hold on;

vals = [];
for I = 1:N,
    valu = x(I);
    height = gammaGuess(I); % generate data set for histogram
    height = round(height);
    for J = 1:height,
        vals = [vals valu];
    end
end

est = gamfit(vals+.00001);
C2 = est(1);
C3 = est(2);

sorig = sum(gammaGuess);
snew  = sum(gampdf(x,C2,C3));
[Yorig,Iorig] = max(gammaGuess);
[Ynew,Inew] = max(gampdf(x,C2,C3));

C1 = 20; %mean(mean(A(3:6,3:6,1)));%sorig/snew;
C4 = Inew-Iorig; %20

%plot(x-30,C1*gampdf(x,C2,C3),'-r');

%% Initialize best guess
Ap = zeros(8,8,4); % pixle loc, pixel loc, params: amplitude, shape1, shape2, delay
Bp = zeros(8,8,N);
for I = 2:7,
    for J = 2:7,
        Ap(I,J,1) = C1; % amplitide
        Ap(I,J,2) = C2; % shape 1
        Ap(I,J,3) = C3; % shape 2
        Ap(I,J,4) = C4; % delay: 0-20
        Bp(I,J,:) = C1*gampdf(x-C4,C2,C3);
    end
end
%figure;imshow(abs(Bp(:,:,60)),[])

%% Generate data based on guess

theta = 0;
for proj = 1:N,
    theta = theta + 111.246;
    im = squeeze(Bp(:,:,proj));
    Rp(:,proj) = radon(im,theta);
end
figure;imshow(Rp,[])

%% Compare errors in the sinogram
error = sqrt(sum(sum((Rp - R).*conj(Rp - R))))
rmsError = dynProjKernel(Ap(:,:,1),Ap(:,:,2),Ap(:,:,3),Ap(:,:,4),R)

figure;imshow(abs([Rp ; R]),[]);

%% modify guess to more closely match "true" sinogram


