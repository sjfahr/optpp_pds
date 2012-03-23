function  difference = dynProjKernel( x, data, npixel)
% input model parameters x create simulated forward projection
% data to compare to the exact data
%  return the vector of the differences between the predicted and measured values

%reshape the model parameters into intuitive matrices
ModelParameters = reshape(x,npixel, npixel,4);
Amplitude = ModelParameters(:,:,1);
Shape     = ModelParameters(:,:,2);
Scale     = ModelParameters(:,:,3);
Delay     = ModelParameters(:,:,4);


N = size(data,2);
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
    Predicted(:,proj) = radon(im,theta); % calculate synthetic acquired sinogram
end

difference = data(:) - Predicted(:);
 
disp(sprintf('rms error %f ',sum(difference.^2)));

end

