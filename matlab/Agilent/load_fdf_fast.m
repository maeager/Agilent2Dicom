function IM = load_fdf_fast(name)
%----------------------------------------
% modified from load_fdf()
%----------------------------------------
% Maj Hedehus, Varian, Inc., Oct 2001.
%----------------------------------------


fid = fopen(name,'r','ieee-le');   % ieee-be for Sun, big endian; ieee-le for Linux, little endian.
if fid == -1
    error(['Can not open file', name]);
end

str = fgetl(fid);
while isempty(findstr('checksum',str))
    if ~isempty(findstr(str,'matrix'))
        if length(findstr(str,','))==2
            V = sscanf(str(findstr(str,'{')+1:length(str)),'%d,%d,%d');
        else
            V = sscanf(str(findstr(str,'{')+1:length(str)),'%d,%d');
        end
    end
    str = fgetl(fid);
end

% Skip past NULL character that separates header and data
v = fread(fid,1,'uchar');
while v ~= 0
    v = fread(fid,1,'uchar');
end

IM = fread(fid,'float');
fclose(fid);
disp 'Reshaping image'
V
IM=reshape(IM,V');
if length(V)==2
    IM=IM';
end
