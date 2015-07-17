function [im,hdr] = readfdf(dirname)
% [im,hdr] = readfdf(dirname)
% dirname is the *.img directory containing the fdf files
% im takes the form [:,:,slice,echo,image]
% hdr contains various info about the data (eg voxelsize) in mm

    files = dir(dirname);
        
    %files = regexp(strtrim(ls('-1', [dirname '/slice*.fdf'])), '\n', 'split');
    idx = cellfun(@(x) ~isempty(x), regexp({files.name}, 'slice'));
    if any(idx)
    
        for f = {files(idx).name}
            f = regexp(f{1}, '/', 'split');
            f = f{end};

            a = sscanf(f, 'slice%03dimage%03decho%03d.fdf');
            slice = a(1);
            image = a(2);
            echo = a(3);

            %[im(:,:,slice, echo, image) hdr.voxelsize hdr.orientation hdr.location] = load_fdf_fast([dirname '/' f]);
            
            [im(:,:,slice,echo,image),fdfhdr] = readfdffile([dirname '/' f]);
            
            if slice == 1
                hdr.voxelsize = fdfhdr.voxelsize;
                hdr.orientation = fdfhdr.orientation;
                hdr.location = fdfhdr.location;
                hdr.TR = fdfhdr.TR;
            end
            
            hdr.TE(echo) = fdfhdr.TE;
            
        end
        return
    end
    
    %files = regexp(strtrim(ls('-1', [dirname '/img_slab*.fdf'])), '\n', 'split');
    idx = cellfun(@(x) ~isempty(x), regexp({files.name}, 'slab'));
    if any(idx)
    
        for f = {files(idx).name}
            f = regexp(f{1}, '/', 'split');
            f = f{end};

            if f(1) == 's'
                a = sscanf(f, 'slab%03dimage%03decho%03d.fdf');
            else
                a = sscanf(f, 'img_slab%03dimage%03decho%03d.fdf');
            end
            slab = a(1);
            image = a(2);
            echo = a(3);

            [im(:,:,:,slab, echo, image) hdr.voxelsize hdr.orientation hdr.location] = load_fdf_fast([dirname '/' f]);
            
        end
        return
    end
    
    im = [];
    

    function [IM voxelsize orientation location] = load_fdf_fast(name)
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
            elseif ~isempty(findstr(str,'roi'))
                ROI = sscanf(str(findstr(str,'{')+1:length(str)),'%f,%f,%f') * 10; %ROI is in cm
            elseif ~isempty(findstr(str,'orientation'))
                orientation = reshape(sscanf(str(findstr(str,'{')+1:length(str)),'%f,%f,%f,%f,%f,%f,%f,%f,%f'),[3 3]);
            elseif ~isempty(findstr(str,'location'))
                location = sscanf(str(findstr(str,'{')+1:length(str)),'%f,%f,%f');
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
        V
        IM=reshape(IM,V');
        if length(V)==2
            IM=IM';
            voxelsize = [ROI(1:2)'./V' ROI(3)];
        else
            voxelsize = ROI'./V';
        end
        
    end
end