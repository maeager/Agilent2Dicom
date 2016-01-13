function [img, hdr, ksp, RE, IM] = readfid(folder)
% [img, hdr, ksp] = readfid(folder)
% read kspace data from Agilent fid and procpar files.
% output img has dimensions [freq, phase, slice, channel, echo]
% note that the image may have to be circular shifted in the phase
% direction.

error(nargchk(1,2,nargin))
 
warning off MATLAB:divideByZero;

% get acqcycles and TE from procpar
[pp hdr] = readpp([folder '/procpar'], true);
% fid = fopen([folder '/procpar'],'r');
% %hdr.acqcycles = []; 
% %hdr.nEchoes = 1;
% while ~feof(fid)
%     line = fgetl(fid);
%     [par,line] = strtok(line);
%     if strcmp(par,'acqcycles')
%         line = fgetl(fid);
%         [~,hdr.acqcycles] = strtok(line);
%         hdr.acqcycles = str2num(hdr.acqcycles);
%     elseif strcmp(par,'TE')
%         line = fgetl(fid);
%         hdr.nEchoes = str2num(strtok(line));
%     end
%     if ~isempty(hdr.acqcycles) && ~isempty(hdr.nEchoes)
%         break
%     end
% end

%fclose(fid);

% open fid file
[fid, errmsg]= fopen(fullfile(folder, 'fid'),'r', 'b');
if fid
[fname,permission,machinefmt,encodingOut] = fopen(fid);

% Read datafileheader
hdr.nblocks   = fread(fid,1,'int32');
hdr.ntraces   = fread(fid,1,'int32');
hdr.np        = fread(fid,1,'int32');
hdr.ebytes    = fread(fid,1,'int32');
hdr.tbytes    = fread(fid,1,'int32');
hdr.bbytes    = fread(fid,1,'int32');
hdr.vers_id   = fread(fid,1,'int16');
hdr.status    = fread(fid,1,'int16');
hdr.nbheaders = fread(fid,1,'int32');

hdr.s_data    = bitget(hdr.status,1,'int16');
hdr.s_spec    = bitget(hdr.status,2,'int16');
hdr.s_32      = bitget(hdr.status,3,'int16');
hdr.s_float   = bitget(hdr.status,4,'int16');
hdr.s_complex = bitget(hdr.status,5,'int16');
hdr.s_hyper   = bitget(hdr.status,6,'int16');



% validate dimensions
if pp.nD == 2
    dims(1) = pp.np/2; % # phase encode lines / 2
    dims(2) = pp.nf/pp.ns; % # frequency lines acquired / # echoes
    dims(3) = pp.ns; % if 2D, # slices, else ni2
    if pp.ni2 > 1 % fse3d sequence has nD == 2, but is a 3d acquisition???
        dims(3) = pp.ni2;
    end
    if regexp(pp.diff{1},'y') > 0
        dims(1) = pp.nphase *(pp.nbdirs+1); % # phase encode lines / 2
        dims(2) = pp.fn1; % # frequency lines acquired / # echoes
        dims(3) = pp.ns; % if 2D, # slices, else ni2
    hdr.np = pp.np;
 %   hdr.s_float=1;
 %   hdr.s_32=1;
    
    end
    
elseif pp.nD == 3
    dims(1) = pp.np/2; % # phase encode lines / 2
    dims(2) = pp.nf/pp.ne; % # frequency lines acquired / # echoes
    dims(3) = pp.ni2; % if 2D, # slices, else ni2  
     if regexp(pp.diff{1},'y')>0
        dims(1) = pp.nphase *(pp.nbdirs+1); % # phase encode lines / 2
        dims(2) = pp.fn1; % # frequency lines acquired / # echoes
        dims(3) = pp.ns; % if 2D, # slices, else ni2
    
    end
else
    error('Can only handle 2D or 3D files (based on procpar field nD)')
end
    
if hdr.np ~= pp.np || ...
   hdr.ntraces ~= pp.nf || ...
   hdr.nblocks ~= pp.arraydim
    error('Cannot resolve fid header with procpar. We''re probably not interpreting the procpar correctly.')
end

% reset output structures
RE = zeros(hdr.np/2, hdr.ntraces, hdr.nblocks, 'single');
IM = zeros(hdr.np/2, hdr.ntraces, hdr.nblocks, 'single');

%% We have to read data every time in order to increment file pointer
nchar = 0;
for i=1:hdr.nblocks  
    fprintf(1, repmat('\b',1,nchar));
    nchar = fprintf(1, 'reading %d values, in block %d of %d',hdr.np*hdr.ntraces, i, hdr.nblocks);
    
    % Read a block header
    scale     = fread(fid,1,'int16');
    bstatus   = fread(fid,1,'int16');
    index     = fread(fid,1,'int16');
    mode      = fread(fid,1,'int16');
    ctcount   = fread(fid,1,'int32');
    lpval     = fread(fid,1,'float32');
    rpval     = fread(fid,1,'float32');
    lvl       = fread(fid,1,'float32');
    tlt       = fread(fid,1,'float32');
    
    hdr.b_data    = bitget(bstatus,1,'int16');
    hdr.b_spec    = bitget(bstatus,2,'int16');
    hdr.b_32      = bitget(bstatus,3,'int16');
    hdr.b_float   = bitget(bstatus,4,'int16');
    hdr.b_complex = bitget(bstatus,5,'int16');
    hdr.b_hyper   = bitget(bstatus,6,'int16');
    hdr.b_more    = bitget(bstatus,8,'int16');
    hdr.b_npcplx    = bitget(bstatus,9,'int16');
    hdr.b_nfcplx      = bitget(bstatus,10,'int16');
    hdr.b_nicplx   = bitget(bstatus,11,'int16');
    hdr.b_ni2cplx = bitget(bstatus,12,'int16');
    hdr.mode_np_ph    = bitget(mode,1,'int16');
    hdr.mode_np_av    = bitget(mode,2,'int16');
    hdr.mode_np_pwr      = bitget(mode,3,'int16');
    hdr.mode_nf_ph    = bitget(mode,5,'int16');
    hdr.mode_nf_av    = bitget(mode,6,'int16');
    hdr.mode_nf_pwr      = bitget(mode,7,'int16');   
    hdr.mode_ni_ph    = bitget(mode,9,'int16');
    hdr.mode_ni_av    = bitget(mode,10,'int16');
    hdr.mode_ni_pwr      = bitget(mode,11,'int16');
    hdr.mode_ni2_ph    = bitget(mode,13,'int16');
    hdr.mode_ni2_av    = bitget(mode,14,'int16');
    hdr.mode_ni2_pwr      = bitget(mode,15,'int16');
    
    
    
    if hdr.s_float == 1
        data = fread(fid,hdr.np*hdr.ntraces,'*float32');
        str='reading floats';
    elseif hdr.s_32 == 1
        data = fread(fid,hdr.np*hdr.ntraces,'*int32');
        str='reading 32bit';
    else
        data = fread(fid,hdr.np*hdr.ntraces,'*int16');
        str='reading 16bit';
    end

    data = reshape(data,hdr.np,hdr.ntraces);
    RE(:,:,i) = data(1:2:hdr.np,:); 
    IM(:,:,i) = data(2:2:hdr.np,:); 
end

fclose(fid);

fprintf(1, '\nPerforming fourier transform...');

%hdr.nChannels = hdr.nblocks/pp.acqcycles;
hdr.nPhaseEncodes = hdr.ntraces/hdr.nEchoes;

%img = zeros(size(RE,1), size(RE,2)/hdr.nEchoes, size(RE,3)/hdr.nChannels, hdr.nChannels, hdr.nEchoes, 'single');
ksp = zeros(dims(1), dims(2), dims(3), hdr.nChannels, hdr.nEchoes, 'single');
img = zeros(dims(1), dims(2), dims(3), hdr.nChannels, hdr.nEchoes, 'single');
kspace = zeros(dims(1), dims(2), dims(3), hdr.nChannels, hdr.nEchoes, 'single');
% img = complex(RE, IM);
% return

if pp.nD == 2 && pp.ni2 == 1
    for echo = 1:hdr.nEchoes
        for channel = 1:hdr.nChannels
            for slice = 1:dims(3)
                %ksp(:,:,slice,channel,echo) =complex(RE(:,slice:dims(3):end,channel), IM(:,channel:dims(3):end,channel));
                ksp(:,:,slice,channel,echo) =complex(RE(:,(slice-1)*dims(2)+1:(slice)*dims(2),channel), IM(:,(slice-1)*dims(2)+1:(slice)*dims(2),channel));
                %ksp(:,:,slice,channel,echo) =complex(RE(:,slice:dims(2)+1:(slice)*dims(2),channel), IM(:,(slice-1)*dims(2)+1:(slice)*dims(2),channel));
                if isfield(pp, 'pelist') && (length(pp.pelist)==size(ksp,2))
                    kspace(:,:,slice,channel,echo)=ksp(:,pp.pelist-min(pp.pelist)+1,slice,channel,echo);
                %ksp(:,:,slice,channel,echo) =complex(RE(:,echo:hdr.nEchoes:end,channel:hdr.nChannels:end), IM(:,echo:hdr.nEchoes:end,channel:hdr.nChannels:end));
                    img(:,:,slice,channel,echo) = fftshift(ifftn(ifftshift(ksp(:,pp.pelist-min(pp.pelist)+1,slice,channel,echo))));
                else
                    kspace(:,:,slice,channel,echo)=ksp(:,:,slice,channel,echo);
                %ksp(:,:,slice,channel,echo) =complex(RE(:,echo:hdr.nEchoes:end,channel:hdr.nChannels:end), IM(:,echo:hdr.nEchoes:end,channel:hdr.nChannels:end));
                    img(:,:,slice,channel,echo) = fftshift(ifftn(ifftshift(ksp(:,:,slice,channel,echo))));
         
                end
                
           end
        end
    end
else %if pp.nD == 3
    for echo = 1:hdr.nEchoes
        for n = 1:hdr.nChannels
            ksp(:,:,:,n,echo) = complex(RE(:,echo:hdr.nEchoes:end,n:hdr.nChannels:end), IM(:,echo:hdr.nEchoes:end,n:hdr.nChannels:end));
            if isfield(pp, 'pelist')
                img(:,:,:,n,echo) = fftshift(ifftn(ifftshift(ksp(:,pp.pelist-min(pp.pelist)+1,:,n,echo))));
            else
                img(:,:,:,n,echo) = fftshift(ifftn(ifftshift(ksp(:,:,:,n,echo))));
            end
        end
    end
end

hdr.pp = pp;

fprintf(1, '\n');
else
    display errmsg
end
