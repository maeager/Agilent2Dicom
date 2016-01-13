function [im] = readraw(dirname,matsizes)
% [im,hdr] = readfdf(dirname)
% dirname is the *.img directory containing the fdf files
% im takes the form [:,:,slice,echo,image]
% hdr contains various info about the data (eg voxelsize) in mm

    files = dir(dirname);
        nchar=0;im=zeros(matsizes);
    %files = regexp(strtrim(ls('-1', [dirname '/slice*.fdf'])), '\n', 'split');
    idx = cellfun(@(x) ~isempty(x), regexp({files.name}, 'slice'));
    if any(idx)
    
        for f = {files(idx).name}
            f = regexp(f{1}, '/', 'split');
            f = f{end};

            a = sscanf(f, 'image%03dslice%03decho%03d.raw');
            slice = a(1);
            image = a(2);
            echo = a(3);
            %fprintf(1, repmat('\b',1,nchar));
            nchar = fprintf(1, 'reading %d values, in slice %d of image %d\n', matsizes(1)*matsizes(2), slice,image);
    
            %[im(:,:,slice, echo, image) hdr.voxelsize hdr.orientation hdr.location] = load_fdf_fast([dirname '/' f]);
            %display(size(im))
            [im(:,:,slice,echo,image)] = readrawfile([dirname '/' f],matsizes);
            
             
        end
        fprintf(1,'\n');
        return
    end
    

    im = [];
    
end