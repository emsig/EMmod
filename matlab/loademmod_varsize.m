function [data,xvec,yvec] = loademmod_varsize(filename,xsize,ysize,dx,dy)
% LOADEMMOD_VARSIZE loads data created with emmod and the corresponding
% coordinate vectors
% 
% Input arguments:
% - filename: path and name of the file to be loaded
% - xsize: number of points in x-direction
% - ysize: number of points in y-direction
% - dx: spacing in x-direction
% - dy: spacing in y-direction
% 
% Output argument:
% - data: matrix containing the data
% - xvec: vector containing the coordinates in the x-direction
% - yvec: vector containing the coordinates in the y-direction
% 
% Usage:
% [data,xvec,yvec] = loademmod(filename,xsize,ysize,dx,dy)

fid = fopen(filename,'r');

% Check if the filesize matches the number of elements specified
status = fseek(fid,0,'eof');
filesize = ftell(fid);
nel = round(filesize/8/2);
if nel~=xsize*ysize
    error('xsize and ysize are not specified properly.')
end
status = fseek(fid,0,'bof');

% Load the data
temp = fread(fid,2*xsize*ysize,'float64',0,'ieee-le');
fclose(fid);
temp2 = reshape(temp,[2,xsize,ysize]);
data = complex(squeeze(temp2(1,:,:)),squeeze(temp2(2,:,:)));
xvec = linspace(-xsize/2,xsize/2-1,xsize)*dx;
yvec = linspace(-ysize/2,ysize/2-1,ysize)*dy;