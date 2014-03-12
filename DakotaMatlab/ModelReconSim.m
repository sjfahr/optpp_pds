% The purpose of this routine is to input the measured MR data
% in this case the MR data is simulated
% An initial guess for the model parameters is created from the data

%% initialize
N = 100; % number of projections
%x = linspace(-30,120,N); % not sure if need negative
time = 1:N;


% setup theoretical phantom
npixel = 16 ; 

ExactDataPyruvate         = zeros(npixel,npixel,4); % pixel loc, pixel loc, params: amplitude, shape, scale, delay
ExactDataPyruvate(:,:,3)  = ones(npixel,npixel); % pixel loc, pixel loc, params: amplitude, shape, scale, delay
ExactDataPyruvate(2:npixel-1,2:npixel-1,1) =    10*rand(npixel-2); % amplitide: 0-10
ExactDataPyruvate(2:npixel-1,2:npixel-1,2) = 1.4 + .2*rand(npixel-2); % shape : 1.4-1.6
ExactDataPyruvate(2:npixel-1,2:npixel-1,3) = 6.5 + 2*rand(npixel-2); % scale : 6.5-8.5
ExactDataPyruvate(2:npixel-1,2:npixel-1,4) = 10+3*rand(npixel-2); % delay: 10-13

ExactDataLactate         = ExactDataPyruvate         
ExactDataLactate(2:(npixel-1),2:(npixel-1),3)  = ExactDataLactate(2:(npixel-1),2:(npixel-1),3)  +10 
ExactDataLactate(2:(npixel-1),2:(npixel-1),4)  = ExactDataLactate(2:(npixel-1),2:(npixel-1),4)  +2

%% Generate synthetic data
theta = 0;
%% initialize image of parameters at a time point
B  = zeros(npixel,npixel); 
for proj = 1:N,
    theta = theta + 111.246; % increment based on golden ratio
    B(:,:) = ExactDataPyruvate(:,:,1).*gampdf(time(proj)-ExactDataPyruvate(:,:,4),ExactDataPyruvate(:,:,2),ExactDataPyruvate(:,:,3)); % bring them into time series
    %B(5:6,6:7,:) = 20*B(5:6,6:7,:); % represent tumor
    MRDataPyruvate(:,proj) = radon(B,theta); % calculate synthetic acquired sinogram
    B(:,:) = ExactDataLactate(:,:,1).*gampdf(time(proj)-ExactDataLactate(:,:,4),ExactDataLactate(:,:,2),ExactDataLactate(:,:,3)); % bring them into time series
    MRDataLactate(:,proj) = radon(B,theta); % calculate synthetic acquired sinogram
end

save('simulatedData.mat','MRDataPyruvate','ExactDataPyruvate','MRDataLactate','ExactDataLactate');
%figure;imshow(MRData,[])

%% recon the images using regular backprojection (optional)

% proj = 1:1:N;
% psi = (proj-1)*111.246;
% imr = iradon(R,psi); % produce the image
% %figure;imshow(imr,[])

