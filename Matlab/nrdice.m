function [dices]=nrdice(epath, scannum ,mpath, isoc)
%nrdice takes two temperature maps and computes the Dice Coefficient as a
%function of time. The script assumes that the two files are spatially and
%temporally registered. Unfiltered by default. 

%example: nrdice('/FUS4/data2/nanorods/20101019nanorods/e7605',22,'/data/cjmaclellan/mdacc/nano/modeltemperaturefull',[5]);
% selects the 22nd scan from the first directory as the first temperature
% array. It selects the vtk files from the second directory and displays
% the 5 degree isotherm.
%dices gives a cell array that gives the dice coefficient at each time point


[~, etmap, fileinfo]=imageread4(epath,scannum); %import experimental tmap

%etmap=wiener3D(etmap);
 
    nrow = double(fileinfo.Rows);
    ncol = double(fileinfo.Columns);
   % te = fileinfo.EchoTime;
   % tr=fileinfo.RepetitionTime
    pixsp=fileinfo.PixelSpacing;
    numtimepts=numel(etmap(1,1,:));
    acqd=fileinfo.Private_0019_105a; %import pixel spacing and timestep

%import model data (same size as experimental data)
%mtmap=readVTK2(mpath,numel(etmap(1,1,:))); 
load(mpath);
mtmap=Model;

pixspx=pixsp(1);
xaxis=(1:nrow)*pixspx/10;
pixspy=pixsp(2);
yaxis=(1:ncol)*pixspy/10;
timestep=acqd/numtimepts/1e6;
taxis=(0:(numel(mtmap(1,1,:))-1))*timestep;

%if experimental data is longer, truncate it to match the model 
if numel(etmap(1,1,:)) > numel(mtmap(1,1,:))
    etmap=etmap(:,:,1:numel(mtmap(1,1,:)));
else
end


mtmapb=mtmap;
etmapb=etmap;
otmapb=etmap;
dice=ones(1,numel(mtmap(1,1,:)));
dices=cell(numel(isoc));
for nn=1:numel(isoc);% for each isotemp line find the matrix that defines the overlapping area
    clevel=isoc(nn);
for kk=1:numel(mtmap(1,1,:));
    for ii=1:numel(mtmap(:,1,1));
        for jj=1:numel(mtmap(1,:,1));
            if mtmap(ii,jj,kk)>=clevel;
                mtmapb(ii,jj,kk)=1;
            else
                    mtmapb(ii,jj,kk)=0;
            end
           if etmap(ii,jj,kk)>=clevel;
                etmapb(ii,jj,kk)=1;
            else
                etmapb(ii,jj,kk)=0;            
           end
             if etmapb(ii,jj,kk)+mtmapb(ii,jj,kk)==2;
                 otmapb(ii,jj,kk)=1;
             else
                 otmapb(ii,jj,kk)=0;
             end
        end
    end

marea(kk)=nnz(mtmapb(:,:,kk));
earea(kk)=nnz(etmapb(:,:,kk));
oarea(kk)=nnz(otmapb(:,:,kk));

if oarea(kk)>0;%if there is overlap compute Dice Co
dice(kk)=2*oarea(kk)/(marea(kk)+earea(kk));
else 
    dice(kk)=0;
end
% imagesc(etmapb(:,:,kk))
% hold on
% contour(etmap(:,:,kk),[5 5])
% figure
% imagesc(mtmapb(:,:,kk))
% hold on
% contour(mtmap(:,:,kk),[5 5])
% figure
% imagesc(otmapb(:,:,kk),[5 5])

end
dices{nn}=dice;
plot(taxis,dices{nn}) %plot dice as a function of time
hold all
end

title(['Dice Coefficient of ' num2str(isoc) ' Degree Isotherms'])
xlabel('Time (s)')
ylabel('Dice Coefficient')

figure

tmax=max(max(max(mtmap,[],3)));
[ii,jj]=find(max(mtmap,[],3)==tmax);
time=find(squeeze(mtmap(ii,jj,:))==tmax)-1;

%  time=input('Timepoint of interest in seconds? '); %ask for specific timepoint
%  time=time/timestep;
 
 imagesc(mtmap(:,:,round(time)));
 %imagesc(mtmapb(:,:,round(time)));
 
 colormap(hot)
 axis off
 axis image
 colorbar
 hold on

 [ccm,hhm] = contour(mtmap(:,:,round(time)),[isoc(1) isoc(1)]);
 [cce,hhe] = contour(etmap(:,:,round(time)),[isoc(1) isoc(1)]);
 
 marea2=marea(round(time))*pixspx*pixspy
 earea2=earea(round(time))*pixspx*pixspy
 

 set(hhm,'color', 'g', 'LineWidth', 3);
 set(hhe,'color', 'b', 'LineWidth', 3);

 title([num2str(isoc) ' Degree Isocontour Lines at ' num2str(round(time*5)) ' Seconds'] );
 legend([hhe,hhm],'Experiment','Model' ); %plot Dice
 
hold off

                
            