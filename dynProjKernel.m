function [ rmsError ] = dynProjKernel( Amplitude,Shape,Scale,Delay,MRData)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


N = size(MRData,2);
x = 1:N;
sDim1 = size(Amplitude,1);
sDim2 = size(Amplitude,2);
Am = repmat(Amplitude,[1 1 N]);
Sh = repmat(Shape,[1 1 N]);
Sc = repmat(Scale,[1 1 N]);
De = repmat(Delay,[1 1 N]);
xx = repmat(x',[1 sDim1 sDim2]); xx = permute(xx,[2 3 1]);


B = Am.*gampdf(xx-De,Sh,Sc); % bring them into time series %% SEE IF MATRIX?????
nans = isnan(B);
B(nans) = 0;

theta = 0;
for proj = 1:N,
    theta = theta + 111.246; % increment based on golden ratio
    im = squeeze(B(:,:,proj));
    R(:,proj) = radon(im,theta); % calculate synthetic acquired sinogram
end

rmsError = sqrt(sum(sum((MRData - R).*conj(MRData - R))));

end

