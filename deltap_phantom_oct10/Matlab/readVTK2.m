function array = readVTK2(vtkfile,timepoints)
%for VTK files starting with 0001. e.g.
%VTKfile.0001.vtk,VTKfile.0002.vtk...
%NOT
%VTKfile.0000.vtk,VTKfile.0001.vtk...


array=zeros(256);
for ii = 1:timepoints
    
%fid = fopen(vtkfile,'r','b');
%sprintf('%s.%04d.vtk',vtkfile,ii)


fid = fopen(sprintf('%s.%04d.vtk',vtkfile,ii),'r','b');

if fid == -1
    return
end

fgetl(fid); % # vtk Datafile Version 3.6
fgetl(fid); % comments
fgetl(fid); % BINARY
fgetl(fid); % DATASET STRUCTURED_POINTS

s = fgetl(fid); % DIMENSIONS nx ny nz
sz = sscanf(s, '%*s%d%d%d').';

fgetl(fid); % SPACING
fgetl(fid); % ORIGIN
fgetl(fid); %POINT_DATA

fgetl(fid); % SCALARS name data_type (eg. SCALARS scalars float)
%svstr = sscanf(s, '%s', 1)
%dtstr = sscanf(s, '%*s%*s%*s')
fgetl(fid); % the LOOKUP_TABLE

V = fread(fid,prod(sz),'float');

V = reshape(V,sz);
array(:,:,ii) = V(:,:);

fclose(fid);
end






