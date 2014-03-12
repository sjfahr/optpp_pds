function [] = writeVTK(vol,vtkfile, fileID,varName,header)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: writeVTK(vol,vtkfile)
%
%   vol:     The 3D matrix to be saved to file
%   vtkfile: The output filename (string)
%   notes:   Only writes binary STRUCTURED_POINTS
%  
% Erik Vidholm 2006

% edited by Josh Yung 5/4/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ischar(vtkfile) == 0
    error('VTKFILE must be a string.')
end

% dimensions
volinfo = whos('vol');

sz = [size(vol,1) , size(vol,2), size(vol,3)] ;
X = sz(1);
Y = sz(2);
Z = sz(3);


% pixel spacing
spX = header.PixelSpacing(1);
spY = header.PixelSpacing(2);
spZ = header.SliceThickness;

for ii = 1:sz(3)
    voltemp = vol;

%vtkfile = 'vtk
fid = fopen(sprintf('%s.%04d.vtk',vtkfile,fileID),'w','b');

% write header
fprintf(fid, '%s\n', '# vtk DataFile Version 3.6');
fprintf(fid, '%s\n', 'created by writeVTK (Matlab implementation by Erik Vidholm)');
fprintf(fid, '%s\n', 'BINARY');  
fprintf(fid, '%s\n', 'DATASET STRUCTURED_POINTS');  
fprintf(fid, '%s%d%c%d%c%d\n', 'DIMENSIONS ', X, ' ', Y, ' ', Z);
fprintf(fid, '%s%f%c%f%c%f\n', 'ORIGIN ', 0.0, ' ', 0.0, ' ', 0.0); 
fprintf(fid, '%s%.4f%c%.4f%c%d\n', 'SPACING ', spX, ' ', spY, ' ', spZ); 
fprintf(fid, '%s%d\n', 'POINT_DATA ', X*Y*Z);

tp = volinfo.class;
datatypeTemplate = 'SCALARS %s %s\n';
if( strcmp(tp, 'uint8') > 0 )
  fprintf(fid, datatypeTemplate, varName,'unsigned_char');
elseif( strcmp(tp, 'int8') > 0 )
  fprintf(fid, datatypeTemplate, varName,'char');
elseif( strcmp(tp, 'uint16') > 0 )
  fprintf(fid, datatypeTemplate, varName,'unsigned_short');
elseif( strcmp(tp, 'int16') > 0 )
  fprintf(fid, datatypeTemplate, varName,'short');
elseif( strcmp(tp, 'uint32') > 0 )
  fprintf(fid, datatypeTemplate, varName,'unsigned_int');
elseif( strcmp(tp, 'int32') > 0 )
  fprintf(fid, datatypeTemplate, varName,'int');
elseif( strcmp(tp, 'single') > 0 )
  fprintf(fid, datatypeTemplate, varName,'float');
elseif( strcmp(tp, 'double') > 0 )
  fprintf(fid, datatypeTemplate, varName,'double');
end

fprintf(fid, '%s\n', 'LOOKUP_TABLE default');

% write data as binary
fwrite(fid,voltemp,tp);

% close file
fclose(fid);
end

