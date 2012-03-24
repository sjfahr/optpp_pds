function  difference = dynProjKernel( x, data, npixel,roi)
% input model parameters x create simulated forward 
% projection data to compare to the exact data
%
%     npixel x npixel is the dimension of the full image
%        roi is the subset of image parameters where 
%        we will acquire the data
%
%  return the vector of the differences between the predicted and measured values

%reshape the model parameters into intuitive matrices
ModelParameters = reshape(x,roi(1,2)-roi(1,1)+1, roi(2,2)-roi(2,1)+1,4);
Amplitude = ModelParameters(:,:,1);
Shape     = ModelParameters(:,:,2);
Scale     = ModelParameters(:,:,3);
Delay     = ModelParameters(:,:,4);


%TODO need to update with acquistion time
N = size(data,2);
time = 1:N;

%initialize image
B = zeros(npixel);

theta = 0;
for proj = 1:N,
    theta = theta + 111.246; % increment based on golden ratio
    B(roi(1,1):roi(1,2),roi(2,1):roi(2,2)) = Amplitude.*gampdf(time(proj)-Delay,Shape,Scale); 
    nans = isnan(B);
    B(nans) = 0;
    Predicted(:,proj) = radon(B,theta); % calculate synthetic acquired sinogram
end

difference = data(:) - Predicted(:);
 
disp(sprintf('rms error %f ',sum(difference.^2)));

end

