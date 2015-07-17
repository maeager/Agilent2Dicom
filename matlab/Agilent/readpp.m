function [pp acq] = readpp(filename, asstruct)
% [pp acq] = readpp(filename, asstruct)
% read Agilent procpar file
% set asstruct to true to return a structure rather than an array
% Optional output acq contains the main acquisition parameters users are
% usually interested in.

error(nargchk(1,2,nargin))

fid = fopen(filename, 'r');

n = 1;
while ~feof(fid)
    line = fgetl(fid);
    [pp(n).name,line] = strtok(line);
    [pp(n).subtype,line] = strtok(line); pp(n).subtype = str2num(pp(n).subtype); 
        %0 (undefined), 1 (real), 2 (string), 
        %3 (delay), 4 (flag), 5 (frequency), 6 (pulse), 7 (integer).
    [pp(n).basictype,line] = strtok(line); pp(n).basictype = str2num(pp(n).basictype);
        %0 (undefined), 1 (real), 2 (string).
    [pp(n).maxvalue,line] = strtok(line); %the maximum value that the parameter can
        %contain, or an index to a maximum value in the parameter parmax 
        %(found in /vnmr/conpar). Applies to both string and real types of 
        %parameters.
    [pp(n).minvalue,line] = strtok(line); %the minimum value that the parameter can 
        %contain or an index to a minimum value in the parameter parmin 
        %(found in /vnmr/conpar). Applies to real types of parameters only.
    [pp(n).stepsize,line] = strtok(line); %a real number for the step size in 
        %which parameters can be entered or index to a step size in the 
        %parameter parstep (found in /vnmr/conpar). If stepsize is 0, it is 
        %ignored. Applies to real types only.
    [pp(n).Ggroup,line] = strtok(line); %0 (ALL), 1 (SAMPLE), 2 (ACQUISITION),
        %3 (PROCESSING), 4 (DISPLAY), 5 (SPIN).
    [pp(n).Dgroup,line] = strtok(line); %The specific application determines the 
        %usage of this integer.
    [pp(n).protection,line] = strtok(line); %a 32-bit word made up of the following 
        %bit masks, which are summed to form the full mask:
        % 0  1    Cannot array the parameter
        % 1  2    Cannot change active/not active status
        % 2  4    Cannot change the parameter value
        % 3  8    Causes _parameter macro to be executed (e.g., if parameter is named sw, the macro _sw is executed when sw is changed)
        % 4  16   Avoids automatic redisplay
        % 5  32   Cannot delete parameter
        % 6  64   System parameter for spectrometer or data station
        % 7  128  Cannot copy parameter from tree to tree
        % 8  256  Cannot set array parameter
        % 9  512  Cannot set parameter enumeral values
        % 10 1024 Cannot change the parameter's group
        % 11 2048 Cannot change protection bits
        % 12 4096 Cannot change the display group
        % 13 8192 Take max, min, step from /vnmr/conpar parameters parmax, parmin, parstep.
    [pp(n).active,line] = strtok(line); %0 (not active), 1 (active).
    [pp(n).intptr,line] = strtok(line); %not used (generally set to 64).
    
    line = fgetl(fid);
    [pp(n).numvalues,line] = strtok(line, ' '); pp(n).numvalues = str2num(pp(n).numvalues);
    
    if pp(n).basictype == 1
        pp(n).values = zeros(1, pp(n).numvalues);
        for i = 1:pp(n).numvalues
            [value,line] = strtok(line, ' ');
            pp(n).values(i) = str2num(value);
        end
        pp(n).lastline = fgetl(fid);
    elseif pp(n).basictype == 2
        pp(n).values = cell(1, pp(n).numvalues);
        pp(n).values{1} = strtrim(line);
        for i = 2:pp(n).numvalues
            pp(n).values{i} = strtrim(fgetl(fid));
        end
        pp(n).lastline = fgetl(fid);
    else
        pp(n).lastline = fgetl(fid);
    end
    
    n = n+1;
end

for p = pp
    pps.(p.name) = p.values;
end

if nargout > 1
    warning 'Acquisition fields may not be correct. Unsure of procpar interpretation.'
    
    % Echo time
    try
        acq.TE = pps.TE;
    catch ME
        acq.TE = pps.te;
    end
    
    % FOV
    acq.volumes = pps.lpe3;
    acq.nEchoes = pps.ne;
    rcvrs = regexp(pps.rcvrs, 'y','match');
    acq.nChannels = length(rcvrs{1});
    acq.mode = sprintf('%dD', pps.nD);
    if pps.nD == 2
        acq.FOVcm = [pps.lro pps.lpe];
        acq.dims = [pps.nf/pps.ns pps.np/2 pps.ns];
        acq.voxelmm = [acq.FOVcm./acq.dims(1:2) pps.thk]*10;
    elseif pps.nD == 3
        acq.FOVcm = [pps.lro pps.lpe pps.lpe2];
        acq.dims = [pps.nf pps.np/2 pps.nv2];
        acq.voxelmm = acq.FOVcm./acq.dims*10;
    end
    
end

% CONVERT OUTPUT TO STRUCT
if nargin == 2 && asstruct
    pp = pps;
end



