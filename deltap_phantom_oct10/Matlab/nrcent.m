function [] = nrcent(path,scannum)
%nr cent takes an array of temperature images and finds the center of
%heating and the offset from the physical center of a cylindrical phantom.
%Note:assumes source is directed along y axis.

%Written by Chris MacLellan 4/11

%Example:nrcent('/FUS4/data2/nanorods/20101019nanorods/e7605',22);
    %selects the 22nd scan from the directory and plots the tmap at max t.
    %Define the center of heating with an ROI. A plot of the
    %center of heating as a function of rows is created and the estimated
    %offset is returned in the command window.An additional plot of the
    %cylinder profile is created from the magnitude image to confirm that
    %the two maximum values represent the cylinder size. 

[dat2, tmap, ~]=imageread4(path,scannum); % import tmap and phase arrays

tmap=wiener3D(tmap,[7,7]);
dat2=wiener3D(dat2,[7,7]);

tmax=max(max(max(tmap,[],3)));
[mi,mj]=find(max(tmap,[],3)==tmax);
mk=find(squeeze(tmap(mi,mj,:))==tmax);
maxpix=[mi, mj, mk]; %find the highest valued pixel in space and time

imagesc(tmap(:,:,mk)); % plot the maximum temperature distribution 
roi=getrect; % The center of heating will be measured in this ROI


side=zeros(1,numel(tmap(:,1,1)));
cent=zeros(1,numel(tmap(:,1,1)));

for kk=(mk-3):(mk+3);
    kk;
    for ii=round(roi(2)):round((roi(2)+roi(4)));
        cent(ii)=cent(ii)+sum(tmap(ii,:,kk).*(1:numel(tmap(1,:,1))))./sum(tmap(ii,:,kk));
        side=side+abs(dat2(ii,:,kk));
    end
end
cent=cent/numel((mk-3):(mk+3));

%^Take a weighted average of the heating profile at each row in the ROI
%over 13 timepoints around max heating (semi-arbitrary 6 frames=30s). Also
%define a vector whose two highest values will give the position of the edges of the
%clyinder

left=find(side==max(side(1:numel(side)/2)));
right=find(side==max(side(numel(side)/2:numel(side))));
middle=(right+left)/2
%^ Find edges of cylinder and physical middle

offset=mode(round(cent(round(roi(2)):round((roi(2)+roi(4))))))-middle
%define the offset, the number of pixels between the physical center and the mode
%of the center of heating. Use the mode because the mean will need to be
%rounded anyways and there can be some bias based on distance from center
%of heating

figure(2)
plot(round(cent))
xlabel('Row')
ylabel('Column')
title('Center of Heating')
figure(3)
plot(side)
xlabel('Column')
ylabel('Magnitude (AU)')
title('Cylinder Profile from Magnitude Image')
%plot the position of the center of heating as a function of the row
%plot cylinder profile from magnitude image. Make sure the two max points
%are indicative of the cylinder width.
        
