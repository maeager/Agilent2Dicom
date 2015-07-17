function [img,hdr] = readfdffile(name)
% [img,hdr] = read_fdf_file(name)

    % open file
    fid = fopen(name,'r','ieee-le');   % ieee-be for Sun, big endian; ieee-le for Linux, little endian.
    if fid == -1
        error(['Can not open file', name]);
    end

    % check endianness - assume first line reads
    % #!/usr/local/fdf/startup
    line = fgetl(fid);        
    if ~strcmp(line, '#!/usr/local/fdf/startup')
        fclose(fid);
        fid = fopen(name,'r','ieee-be');
        if fid == -1
            error(['Can not open file', name]);
        end
        line = fgetl(fid);        
        if strcmp(line, '#!/usr/local/fdf/startup')
            fclose(fid);
            error('File ''%s'' not recognised as fdf format', name)
        end
    end

    % read in header
    line = fgetl(fid);
    while isempty(line) || double(line(1)) ~= 12
        if ~isempty(line) && ~strcmp(line(1:2),'/*')
            [~,b] = strtok(line);
            [a,b] = strtok(b,'=');
            a = strtrim(a);
            if strcmp(a(1),'*')
                a = a(2:end);
            end
            if strcmp(a(end-1:end),'[]')
                a = a(1:end-2);
            end

            b = strtrim(b(2:end));
            b = b(1:end-1);
            if any(strfind(b, '"'))
                b = strrep(b,'"', '''');
            else
                b = strrep(b,'{','[');
                b = strrep(b,'}',']');
            end
            eval(sprintf('hdr.%s = %s;',a,b))
        end
        line = fgetl(fid);
    end

    % Skip past NULL character that separates header and data
    v = fread(fid,1,'uchar');
    while v ~= 0
        v = fread(fid,1,'uchar');
    end

    % read in data
    img = fread(fid,[hdr.storage num2str(hdr.bits)]);

    % close file
    fclose(fid);

    % reshape data
    img = reshape(img,hdr.matrix);
    if hdr.rank == 2
        img = img';
        hdr.voxelsize = hdr.roi ./ [hdr.matrix 1];
    else
        hdr.voxelsize = hdr.roi ./ hdr.matrix;
    end

end
