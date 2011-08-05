function [movz,movt,mtemp,etemp]=nrprofiles(epath, scannum ,mpath)
%nrprofiles takes two temporally and spatially registered temperature arrays
%and generates axial/transverse profiles and plots the temperature as a
%funtion of time at a user specified point.

%example:nrprofiles('/FUS4/data2/nanorods/20101019nanorods/e7605',22,'/data/cjmaclellan/mdacc/nano/modeltemperaturefull');
% selects the 22nd scan from the first directory as the first temperature
% array. It selects the vtk files from the second directory as the second
% temperature array.

%etemp=a vector that gives the experimental temperature as a fucntion of
%time
%mtemp= a vector that gives the model temperature as a function of time
%movt/movz=movies showing the the heating in the transverse and axial
%directions

% written by Chris MacLellan 4/11
[~, etmap, fileinfo]=imageread4(epath,scannum);

    nrow = double(fileinfo.Rows);
    ncol = double(fileinfo.Columns);
   % te = fileinfo.EchoTime;
   % tr=fileinfo.RepetitionTime
    pixsp=fileinfo.PixelSpacing;
     numtimepts=numel(etmap(1,1,:));
    acqd=fileinfo.Private_0019_105a;
%     lon=6
%     loff=42
%     int=
    

%import model data (same size as experimental data)
%mtmap=readVTK2(mpath,numel(etmap(1,1,:)));
load(mpath);
mtmap=Model;

    
% for qq=1:numel(mtmap(1,1,:))-2;
% mtmap2(:,:,qq)=mtmap(:,:,2*qq);
% end
% mtmap=mtmap2;


% for q=1:numel(mtmap(1,1,:))-1;
% mtmap2(:,:,q)=(mtmap(:,:,q)+mtmap(:,:,q+1))/2;
% end
%mtmap=mtmap2;


%if experimental data is longer, truncate it to match the model 
if numel(etmap(1,1,:)) > numel(mtmap(1,1,:))
    etmap=etmap(:,:,1:numel(mtmap(1,1,:)));
else
end

imagesc(max(mtmap,[],3))

tmax=max(max(max(mtmap,[],3)));
[ii,jj]=find(max(mtmap,[],3)==tmax);
imnum=find(squeeze(mtmap(ii,jj,:))==tmax)-1;%not at max time t because of errors associated with laser shutoff time

% [jj,ii]=getpts;
% ii=round(ii);
% jj=round(jj);

% point=input('Point?') ;
% ii=point(1);
% jj=point(2);

if scannum==12;
    int=81;
elseif scannum==22;
    int=73;
else 
    int=ii;
end

pixspx=pixsp(1);
xaxis=((1:nrow)-jj)*pixspx;
%xaxis=1:nrow;
pixspy=pixsp(2);
yaxis=((1:ncol)-int)*pixspy;
%yaxis=1:nrow;
timestep=acqd/numtimepts/1e6;
taxis=(0:(numel(mtmap(1,1,:))-1))*timestep;

clf

for nn=1:numel(mtmap(1,1,:))
     plot(yaxis,mtmap(:,jj,nn),'r');
     hold on;
     plot(yaxis,etmap(:,jj,nn),'b');
     hold off;
    % axis([0,numel(mtmap(:,1,10)),-.5,max(max(max(etmap,[],3),[],2),[],1)])
     ylim([-.5,max(max(max(etmap,[],3),[],2),[],1)])
     movz(nn)=getframe;
     clf
     plot(xaxis,mtmap(ii,:,nn),'r');
     hold on;
     plot(xaxis,etmap(ii,:,nn),'b');
     hold off;
     ylim([-.5,max(max(max(etmap,[],3),[],2),[],1)])
     %axis([0,numel(mtmap(1,:,10)),-.5,max(max(max(etmap,[],3),[],2),[],1)])
     movt(nn)=getframe;
     clf
end
 
mtemp=squeeze(mtmap(ii,jj,:));
etemp=squeeze(etmap(ii,jj,:));
%mtemp2=squeeze(mtmap2(ii,jj,:));

%imnum=input('Image number?')
clf
plot(xaxis,mtmap(ii,:,imnum),'r')
hold;
plot(xaxis,etmap(ii,:,imnum),'b+')
plot(xaxis(jj),mtmap(ii,jj,imnum),'go');
%title('Axial Profile')
ylim([0,15]);
xlim([-15,15]);
xlabel('Distance (mm)')
ylabel('Temperature Increase ^o C')
hold off

figure
plot(yaxis,mtmap(:,jj,imnum),'r')
hold;
plot(yaxis,etmap(:,jj,imnum),'b+')
plot(yaxis(ii),mtmap(ii,jj,imnum),'go');
%temp=mean(abs(dat2),3);
%imagesc(temp)
%plot(yaxis,temp(:,jj),'k')
%title('Depth Profile')
xlim([-30,30]);
ylim([0,15]);
xlabel('Distance (mm)')
ylabel('Temperature Increase ^o C')
hold off

figure
plot(taxis,mtemp,'r')
hold;
plot(taxis,etemp,'b+')
plot(taxis(imnum),mtmap(ii,jj,imnum),'go')
%plot(taxis(1:58)+2.5,mtemp2,'go')
%title('Temperature as Function of Time')
ylim([0,15]);
xlim([0,290.3]);
xlabel('Time (s)')
ylabel('Temperature Increase ^o C')
hold off

%whos