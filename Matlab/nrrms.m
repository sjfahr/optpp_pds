function [diff, rms]=nrrms(epath, scannum ,mpath)
% nrrms finds the rms difference between two previously registered tmaps.

%diff is a matrix defining the difference  between the tmaps. A positive
%entry in diff means the model is hotter than experiment.

%rms is the rms difference between the tmaps

% epath is the path of the experimental data (string)
% scannum is the number of the scan of interest in the path epath
% the file defined by epath and scanum is assumed to be a DICOM file

% mpath is the path of the vtk model file (string). 

%for example: nrrms('/FUS4/data2/nanorods/20101019nanorods/e7605',22,'/data/cjmaclellan/mdacc/nano/modeltemperature');
% loads the 22nd scan from the directory ...20101019 nanorods and loads the
% files modeltemperature.XXXX.vtk from .../mdacc/nano
% Chris MacLellan 3/22/11

%import the experimental data using a modified version of Andrews imageread
%code. By default there is no filtering in imageread4.m

[dat2, etmap, fileinfo]=imageread4(epath,scannum);
 
    nrow = double(fileinfo.Rows);
    ncol = double(fileinfo.Columns);
    te = fileinfo.EchoTime;
   % tr=fileinfo.RepetitionTime
   % pixsp=fileinfo.PixelSpacing;
    pixsp=fileinfo.PixelSpacing;
    numtimepts=fileinfo.NumberOfTemporalPositions;
    acqd=fileinfo.Private_0019_105a;
 
    

%import model data (same size as experimental data)
%mtmap=readVTK2(mpath,numel(etmap(1,1,:)));
load(mpath);
mtmap=Model;
    
pixspx=pixsp(1,1)/10;
xaxis=(1:nrow)*pixspx;
pixspy=pixsp(2,1)/10;
yaxis=(1:ncol)*pixspy;
timestep=acqd/numtimepts/1e6;
taxis=(0:(numel(mtmap(1,1,:))-1))*timestep;

%if experimental data is longer, truncate it to match the model 
if numel(etmap(1,1,:)) > numel(mtmap(1,1,:))
    etmap=etmap(:,:,1:numel(mtmap(1,1,:)));
else
end

%find the difference between experiment and model
diff=mtmap-etmap;

rms=sqrt(sum(diff.^2,3)/numel(diff(1,1,:)));
imagesc(xaxis,yaxis,rms)
axis off
%title('RMS Difference ^oC')
colorbar('westoutside');

% rms2=sqrt(sum(sum(diff.^2,2),1)/numel(diff(:,:,1)));
% rms2=squeeze(rms2);
% figure
% plot(taxis,rms2);
% xlabel('Time (s)')
% ylabel('RMS Difference')
% title('RMS Difference in Time')

%to trace the error in time we must throw out pixel that have poor snr

%find where max heating occurs to calculate signal
tmax=max(max(max(etmap,[],3)));
[jj,ii]=find(max(etmap,[],3)==tmax);
kk=find(squeeze(etmap(jj,ii,:))==max(max(max(etmap,[],3))));

%move roi away from interface
jj=jj+10;

%establish ROI to measure signal in
m=5;
roi=[ii-m,jj-m,2*m,2*m];
figure
imagesc(mean(abs(dat2(:,:,2:4)),3))
hold;
scatter([ii-m,ii-m,ii+m,ii+m],[jj-m,jj+m,jj+m,jj-m],'k+')

mag=abs(dat2);
%find avg signal in area near heating
%[~,avg]=meanroi(imcrop2(mean(mag(:,:,2:4),3),roi));

%find roi to measure noise in
shift=round(35/2/pixsp(1));
roi2=[roi(1)+shift,roi(2),roi(3),roi(4)];
scatter([ii-m+shift,ii-m+shift,ii+m+shift,ii+m+shift],[jj-m,jj+m,jj+m,jj-m],'k+')

[std,~]=meanroi(imcrop2(mean(mag(:,:,2:4),3),roi2))
hold off

%error in signal magnitude = std/.655 for rician dist
sigmaa=std/.655;
alpha=.01;
gamma= 42.576; %Mhz/T;
b=1.5;%T;
te=te*10^-3;%s

%eliminate pixels where signal does not exceed 2 sigmaa

etmap2=etmap;
sigmat=etmap;
mtmap2=mtmap;
bin=etmap-etmap;
for pp=1:numel(etmap(1,1,:));
for nn=1:numel(etmap(:,1,1));
    for mm=1:numel(etmap(1,:,1));
         
        sigmat(nn,mm,pp)=(sqrt(2)*sigmaa)/(2*pi*alpha*(gamma)*b*te*(mag(nn,mm,pp)));
         
        if sigmat(nn,mm,pp)==inf;
        sigmat(nn,mm,pp)=0;
        end
        
        if etmap(nn,mm,pp) <= 3*sigmat(nn,mm,pp);
            etmap2(nn,mm,pp)=0;
            mtmap2(nn,mm,pp)=0;    
        else
            bin(nn,mm,pp)=1;
           etmap2(nn,mm,pp)=etmap(nn,mm,pp);
           mtmap2(nn,mm,pp)=mtmap(nn,mm,pp);
%             etmap2(nn,mm,pp)=1;
%            mtmap2(nn,mm,pp)=1;
        end
    end
end
end

% figure
% montage0(etmap2)
% 
% nomask=nnz((sum(abs(etmap),3)))
% for ss=1:5;
%     gps=nnz(etmap2(:,:,ss))/nomask
% end

diff2=mtmap2-etmap2;
deltat2=(bin).*sigmat;
% figure
% montage0(abs(diff2)>3*deltat2);
% montage0(deltat2)
% figure

rms2=squeeze(diff2(1,1,:));
ut=rms2;
for qq=1:numel(diff2(1,1,:))
rms2(qq)=sqrt(sum(sum(diff2(:,:,qq).^2,2),1)/nnz(diff2(:,:,qq)));
ut(qq)=sqrt(sum(sum(deltat2(:,:,qq).^2,2),1)/nnz(deltat2(:,:,qq)));
end

rms2=squeeze(rms2);
ut=squeeze(ut);
figure
plot(taxis(12:59),rms2(12:59));
% hold;
% plot(taxis,ut);
xlabel('Time (s)')
ylabel('RMS Difference')
%title('RMS Difference in Time')

