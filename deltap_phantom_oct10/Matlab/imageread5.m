function [dat2,tempmap,fileinfo]=imageread4(dname,scannum)
%Modified version of imageread2.m used for generating tmaps. No filtering
%by default.
%Modififed by Chris MacLellan 3/22/11
%Imageread.m belongs in the directory above the image directory 
% containing the dated directory (e.g above the directory 070131)
% This program reads in and processes REAL and IMAGINARY files 
% of 2D images into into a temperature map.
 
% 
% Imageread takes 2 arguments:
%               1) The directory path form the current directory 
%                  level down to but not including the subdirectory
%                  that contains all the scans.
%               2) The second argument is the scan number.

% For example to read the temperature data that is stored in the directory 
% '/FUS4/data2/mfgre/070618_canine_liver/e467461/s468041' on venus you
% would type imageread('/FUS4/data2/mfgre/070618_canine_liver/e467461/',43)
% because s468041 is the 43rd file in the diectory 
% '/FUS4/data2/mfgre/070618_canine_liver/e467461/s468041'.

% If you want to change this please feel free to do so as I tend to write
% crazy programs that not even I can understand.

% Written by Andrew Elliott around Janurary 2007.

path(path,'/FUS/matlab/')


s = dir(dname);                     % Returns information on the files is the directory
                                    % as a structure. 
scan_dir = s(scannum+2).name;
full_dir = [dname,'/',scan_dir,'/'];
s2 = dir(full_dir);
filename = s2(3).name;              % Selects third file in the directory, the first
                                    % and second files are just the '.' and
                                    % '..' files.

% _____________________________________________________________
% Beginning of Section 1
% This point in the script begins reading in for DICOM files.
% For older GE style files this section will have to be rewritten to
% accomedate the readgelx() function.

    if 1
        path(path,full_dir);
    fileinfo = dicominfo(filename);
    %fileinfo = dicominfo(filename,'dictionary','gems-dicom-dictHDX.txt');
    nrow = fileinfo.Rows;
    ncol = fileinfo.Columns;
    te = fileinfo.EchoTime;
    numoffile = fileinfo.ImagesInAcquisition;
    alpha = -.01%-0.0097;
    firstnum = str2num(filename(1,2:regexp(filename,'.MRDC.')-1));
    lastnum = 1;
    ii = 1;
    rawimage = zeros(nrow,ncol,numoffile);
    tempim = zeros(nrow,ncol,2);
    tmap = zeros(nrow,ncol,numoffile/2);
    dat = zeros(nrow,ncol,numoffile/2);
    dat2 = zeros(nrow,ncol,numoffile/2);
    timage = tmap;
    tempmap = tmap;
    
    
   slicenum = 0;    % 2D imaging
  

    for tt = 0:2:numoffile-1
        filename1 = ['i' num2str(firstnum+tt+slicenum) '.MRDC.' num2str(lastnum+tt+slicenum)];
        filename2 = ['i' num2str(firstnum+tt+1+slicenum) '.MRDC.' num2str(lastnum+tt+1+slicenum)];
        tempim(:,:,1) = dicomread(filename1);
        tempim(:,:,2) = dicomread(filename2);
        rawimage(:,:,ii) = cplximg(tempim);
        ii = ii+1;
  %      slicenum = 4*2-1;
    end
    end
    % This ends the part of the script that reads DICOM files.  From this point
    % on all data is in complex form in the 'rawimage' array.
    % End of Section 1
    % ____________________________________________________________

    if 1
           
        % Beginning of Section 2
        % Puts 'rawimage' array into 'dat' array. This section is not necessary, I
        % put it in if you need to manipulate only a subset of the images  
        jj = 1;
        image_num = 2;
        for ii = image_num:numoffile/2
            dat(:,:,jj) = rawimage(:,:,ii); %wiener2(rawimage(:,:,ii),[11,11]);
            jj = jj+1;
        end

        
        % End of Section 2.
        
        % Section 3
        % This section creates an image mask.  The mask eliminates the noise 
        % from the area surrounding the image, making it easier to see.  This
        % mask is used in section 5 to clean up the image.
%          scale_fac = 0.30;
%          msk = mean(rawimage(:,:,1:2),3);
%          mask_level = scale_fac*max(max(abs(msk)));
%          msk = (abs(msk)>mask_level);
%          
%          for tt = 1:(numoffile/2)
%          dat2(:,:,tt) = dat(:,:,tt);
%          end

        % End of Section 3

        % Beginning of Section 4
        % This is where the temperature map is calculated form the complex
        % array 'dat2'.

        for jj = 1:(numoffile/2)-1
        tmap(:,:,jj+1) = ((tmap(:,:,jj)+phitmap(angle(dat2(:,:,jj).*conj(dat2(:,:,jj+1))),te,alpha)));
        end

        % End of section 4
        % The temperature map tmap exists and can be saved
        % at this point.
    end

        % Beginning of Section 5
        % This section cleans up the temperature map and makes it look nice.
        % It does this by first applying a wiener filter to the tmap data. Then
        % it applys the mask 'msk' created in section 4.  After this the
        % maximum temperature image is displayed.
    if 1
        for jj = 1:(numoffile/2)
            %timage(:,:,jj) = wiener2(tmap(:,:,jj).*msk,[3,3]);
            %tempmap(:,:,jj) = wiener2(tmap(:,:,jj).*msk,[5,5]);
            tempmap(:,:,jj) = tmap(:,:,jj);
        end
        
%         a = 1;
%         b = [1/4 1/4 1/4 1/4];
%         for jj = 1:256
%            for ii = 1:256
%                tempmap(ii,jj,:) = filter(b,a,timage(ii,jj,:));
%            end
%         end

        %h1 = figure(1)
        %imagesc(max(tempmap,[],3))
        %colormap 'hot'
        %caxis([0 20])
    end
    
   % roi1 = round(getrect());
   % ttime = tempmap(roi1(1,2):roi1(1,2)+roi1(1,4),roi1(1,1):roi1(1,1)+roi1(1,3),:);
   % y = squeeze(mean(mean(ttime,1),2));
   % err = squeeze(std(std(ttime,[],1),[],2));
%         for jj = 1:(numoffile/2)-1
%             timage(:,:,jj) = (tmap(:,:,jj).*msk);
%         end
%             % Compute a magnitude image
%              h2 = figure(2)
%              magimage = abs(dat2(:,:,1)+dat2(:,:,2)+dat2(:,:,3)+dat2(:,:,4))/4;
%              imagesc(magimage)
%              colormap 'gray'
%         % End of Section 5
   

end

