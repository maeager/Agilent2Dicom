function [kspace,np,ntraces,nblocks] = readVarian2D(filename)

   file_id=fopen(filename, 'r', 'b');		%always big endian format, even for Linux

   %Start reading the header to get the first 6 elements
   [data_holder, count] = fread(file_id, 6, 'int32');
   
   nblocks = data_holder(1); %data blocks in a file (second dimension)	
   
   ntraces = data_holder(2); %slices (can be third direction)	
   
   np = data_holder(3);   % real and imaginary points (first dimension)
   
   
   ebytes = data_holder(4); % 4 bytes per point	
   tbytes = data_holder(5);  %np*ebytes
   bbytes = data_holder(6);  %ntraces*tbytes+nhead*sizeof(datablockhead=28bytes)-total
   
   %Read to see the data type and the status of the binary
   
   [data_holder, count] = fread(file_id, 2, 'int16');
   vers_file_id = data_holder(1);
   status = data_holder(2);
   [data_holder, count] = fread(file_id, 1, 'int32');
   nhead = data_holder(1);
   
   
 
   data_holder = zeros(ntraces*np, 1);
						%determine the binary format from the status info
   BinStat = fliplr(dec2bin(status));
   
   %This is a case of int16 data type
   if ((BinStat(3) == '0') & (BinStat(4) == '0'))
     for nblocks_count=1:nblocks
       for nhead_count=1:nhead
         fread(file_id, 14, 'int16');			%read in the block headers (28 bytes)
       end 
       [b, count] = fread(file_id, ntraces*np, 'int16');   	%read in the actual data (ntraces*np)
       							%put the data in the vector
       data_holder( ((nblocks_count-1)*np*ntraces+1) : nblocks_count*np*ntraces ) = b;
     end  
     
     
   %This is a case of int32 data type
   elseif ((BinStat(3) == '1') & (BinStat(4) == '0'))
     for nblocks_count=1:nblocks
       for nhead_count=1:nhead
          [d,count]=fread(file_id, 14, 'int16');			%read in the block headers (28 bytes)
          bhead=d;
       end 
       [b, count] = fread(file_id, ntraces*np, 'int32');   	%read in the actual data (ntraces*np)
       							%put the data in the vector
       data_holder( ((nblocks_count-1)*np*ntraces+1) : nblocks_count*np*ntraces ) = b;
     end  
     
     
  %This is a case of float data type   
   else  
     for nblocks_count=1:nblocks
       for nhead_count=1:nhead
         fread(file_id, 14, 'int16');			%read in the block headers (28 bytes)
       end 
       [b, count] = fread(file_id, ntraces*np, 'float');   	%read in the actual data (ntraces*np)
       							%put the data in the vector
       data_holder( ((nblocks_count-1)*np*ntraces+1) : nblocks_count*np*ntraces ) = b;
     end  
  end
  size3=ntraces;
  size1=nblocks;
  size2=np;
  
  kspace=complex(zeros(nblocks,ntraces*np/2));
  
  
  for nblocks_count=1:nblocks
     kspace(nblocks_count, 1:(size2/2)*size3) = (data_holder((nblocks_count-1)*size2*size3 + 1 : 2 : nblocks_count*size2*size3-1)+....
     + sqrt(-1)*data_holder((nblocks_count-1)*size2*size3 + 2 : 2 : nblocks_count*size2*size3))';
  end; 
 
  
  
 
