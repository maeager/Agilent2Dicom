function [Rimg, Iimg] = agilentfid(filename)
% m-file that can open Varian FID imaging files in Matlab.
% Usage: img = fdf;
% Your image data will be loaded into Rimg and Iimg
%

warning off MATLAB:divideByZero;

[fid] = fopen(filename,'r', 'b'); %[pathname filename],'r');

% skip file header
fseek(fid, 32, 'bof');

sz = [512 512 512];
Rimg = single(zeros(sz));
Iimg = single(zeros(sz));

for n = 1:sz(3)
    % skip block header 
    fseek(fid, 28, 'cof');
    
    % read in block
    tmp = fread(fid, sz(1)*sz(2)*2, '*float32');
    
    % split interleave
    Rimg(:,:,n) = reshape(tmp(1:2:end),sz(1:2));
    Iimg(:,:,n) = reshape(tmp(2:2:end),sz(1:2));
end

fclose(fid);
