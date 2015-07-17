function [final_kspace,np,ntraces,nblocks,ns]= ReadVarStart
%This file reads in Varian raw fid file
%Output:final_kspace,np,ntraces,nblocks
%NB: this is only for seqcon='ncsnn', 'nccnn' (2D) and 'nccsn' (3D)
%Example:
%[final_kspace,np,ntraces,nblocks]= ReadVarStart(3);

global pathname np ntraces nblocks array seqcon 
%global kspace
global image_space

if ischar(pathname)
    [filename,pathname]=uigetfile('*.*','Select fid file for 2D or 3D',pathname);
else
[filename,pathname]=uigetfile('*.*','Select fid file for the 2D or 3D');
end
filename
if filename ~=0
    %Open the input file
    filename=[pathname filename];
    [fid,msg]=fopen(filename,'rt');
       if fid<0
        % There is an error
        str='File Open Failed'
       else
        filename
       end
else
    error('File Open Stopped')
    final_kspace=0;
    np=0;
    ntraces=0;
    nblocks=0;
    
end

 [kspace,np,ntraces,nblocks] = ReadVarian2D(filename);
 
%seqcon='nccnn' %this is for gems
%seqcon='ncsnn'; %this is for sems
filenameProc=[pathname 'procpar']
   [fid,msg]=fopen(filenameProc,'rt');
    if fid<0
        
     %There is an error
     str=['File Open Failed'];
     errordlg(str,title,'modal');
     
     else
   
  %   sw = ReadProcpar('sw',filenameProc)
   %  sfrq = ReadProcpar('sfrq',filenameProc)
     disp(sprintf('Reading in Varian Procpar files from %s',filename));
     seqcon = ReadProcpar('seqcon',filenameProc)
     ns=ReadProcpar('ns',filenameProc);
     arraydim=ReadProcpar('arraydim',filenameProc);
     
    end       
 
 
if seqcon=='ncsnn'
    
 final_kspace=reshape(kspace,[nblocks np/2 ntraces]);
 
 image_space=abs(fftshift(fftn(final_kspace)));
 %this is an example how to view image_space
 %figure;pcolor(squeeze(image_space(:,:,1)));shading interp;colormap(gray);
 
%  figure;pcolor(abs(fftshift(fft2(squeeze(final_kspace(:,:,nslice))))));
elseif seqcon =='nccnn'
     %disp('need to know the number of slices from the procpar, ns')
     %ns=7; %enter it here if you are not getting it from the procpar
     
     if arraydim>1 %array goes into nblocks
         final_kspace=reshape(kspace,[nblocks np/2 ntraces]);
%          figure;pcolor(abs(fftshift(fft2(squeeze(final_kspace(1,:,:))))));
     else
         
         final_kspace=reshape(kspace,[np/2 ns ntraces/ns]);
%           figure;pcolor(abs(fftshift(fft2(squeeze(final_kspace(:,ns,:))))));
     end
     
     image_space=abs(fftshift(fftn(final_kspace)));
     
elseif seqcon=='nccsn'
    disp('This is 3D case')
    global pelist petable
    petable=ReadProcpar('petable',filenameProc)
    etl=ReadProcpar('etl',filenameProc)
    pelist=ReadProcpar('pelist',filenameProc)
     
    %Check to see if petable is required to reconstruct the data
    if ~isempty(petable)&& (~isempty(etl))>1
      disp('This requires petable in the same folder as the fid data, if pelist doesnt exist')
          if size(pelist,2)>1
          table=pelist.';
          else
      
          table=textread([pathname petable],'%d','delimiter',' ','headerlines',1);
          end
      rawdata=zeros(nblocks,ntraces,np/2);
      
       sstart=1;
 
       for ivv=1:ntraces/etl
          for iecho=1:etl
                 view=(ivv-1)*etl + iecho;
                 iv=table(view)+ntraces/2; %This is because table goes from -pe/2 to +pe/2
                 final_kspace(1:nblocks,iv,:)=kspace(1:nblocks,sstart:sstart+np/2-1);
                 sstart=sstart+np/2;
          end
       end
    else
    
    
    
      final_kspace=reshape(kspace,[nblocks np/2 ntraces]);
    
    end
    image_space=abs(fftshift(fftn(final_kspace)));
    
%     figure;pcolor(squeeze(image_space(nblocks/2,:,:)));shading interp;colormap(gray);
%     figure;pcolor(squeeze(image_space(:,np/4,:)));shading interp;colormap(gray);
%     figure;pcolor(squeeze(image_space(:,:,ntraces/2)));
else
    disp('Not implemented yet')
    return
end
 shading interp;colormap(gray);