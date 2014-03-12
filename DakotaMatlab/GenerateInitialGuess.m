function  InitialGuess = GenerateInitialGuess( data, npixel, roi)

%% make a guess to starting gamma
gammaGuess = sum(data,1);
%figure;plot(x,gammaGuess);
%hold on;

N = size(data,2);
time = 1:N;

vals = [];
for I = 1:N,
    valu = time(I);
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
snew  = sum(gampdf(time,C2,C3));
[Yorig,Iorig] = max(gammaGuess);
[Ynew,Inew] = max(gampdf(time,C2,C3));

C1 = 20; %mean(mean(ExactData(3:6,3:6,1)));%sorig/snew;
C4 = Inew-Iorig; %20

%plot(x-30,C1*gampdf(x,C2,C3),'-r');

%% Initialize best guess
FullImage = zeros(npixel,npixel,4); % pixle loc, pixel loc, params: amplitude, shape1, shape2, delay
FullImage(:,:,1) = C1; % amplitide
FullImage(:,:,2) = C2; % shape 
FullImage(:,:,3) = C3; % scale
FullImage(:,:,4) = C4; % delay: 0-20
%figure;imshow(abs(Bp(:,:,60)),[])

InitialGuess  = FullImage(roi(1,1):roi(1,2),roi(2,1):roi(2,2),:);
end 
