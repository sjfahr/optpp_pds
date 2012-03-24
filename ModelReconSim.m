% The purpose of this routine is to input the measured MR data
% in this case the MR data is simulated
% An initial guess for the model parameters is created from the data

%% initialize
N = 100; % number of projections
%x = linspace(-30,120,N); % not sure if need negative
time = 1:N;


% setup theoretical phantom
npixel = 8 ; 

ExactData         = zeros(npixel,npixel,4); % pixel loc, pixel loc, params: amplitude, shape, scale, delay
ExactData(:,:,3)  = ones(npixel,npixel); % pixel loc, pixel loc, params: amplitude, shape, scale, delay
ExactData(2:npixel-1,2:npixel-1,1) =    10*rand(npixel-2); % amplitide: 0-10
ExactData(2:npixel-1,2:npixel-1,2) = 2 + 5*rand(npixel-2); % shape : 2-7
ExactData(2:npixel-1,2:npixel-1,3) = 1 + 2*rand(npixel-2); % scale : 1-3
ExactData(2:npixel-1,2:npixel-1,4) =    40*rand(npixel-2); % delay: 0-40
ExactData(5:6,6:7,1) = 20*ExactData(5:6,6:7,1); % change model amplitide to represent tumor

%% Generate synthetic data
theta = 0;
%% initialize image of parameters at a time point
B  = zeros(npixel,npixel); 
for proj = 1:N,
    theta = theta + 111.246; % increment based on golden ratio
    B(:,:) = ExactData(:,:,1).*gampdf(time(proj)-ExactData(:,:,4),ExactData(:,:,2),ExactData(:,:,3)); % bring them into time series
    %B(5:6,6:7,:) = 20*B(5:6,6:7,:); % represent tumor
    MRData(:,proj) = radon(B,theta); % calculate synthetic acquired sinogram
end

save('simulatedData.mat','MRData','ExactData');
%figure;imshow(MRData,[])

%% recon the images using regular backprojection (optional)

% proj = 1:1:N;
% psi = (proj-1)*111.246;
% imr = iradon(R,psi); % produce the image
% %figure;imshow(imr,[])

