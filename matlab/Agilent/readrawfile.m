function [img] = readrawfile(name,matsizes)
% [img,hdr] = read_fdf_file(name)

    % open file
    fid = fopen(name,'r');   % ieee-be for Sun, big endian; ieee-le for Linux, little endian.
    if fid == -1
        error(['Can not open file', name]);
    end

   
    % read in data
    img = fread(fid,matsizes(1)*matsizes(2)*2,'*float32');

    % close file
    fclose(fid);

    % reshape data
    img = reshape(img,matsizes);


end
