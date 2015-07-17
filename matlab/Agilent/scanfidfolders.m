function out = scanfidfolders(parent)

    folders = dir([parent '/*.fid']);
    
    datasets = [];
    
    for folder = folders'
        datasets(end+1).folder = folder.name;
        fid = fopen([folder.name '/text'],'r');
        if fid == -1
            error('Cannot open %s', [folder.name '/text'])
        end
        datasets(end).protocol = fgetl(fid);
        fclose(fid);
    end
    
    format = sprintf('%%-%ds%%-%ds\n',5 + max(cellfun(@length,{datasets.folder})), ...
                                  5 + max(cellfun(@length,{datasets.protocol})));
           
    fprintf(1,['\n\n' format], 'FOLDER', 'PROTOCOL');
    for d = datasets
        fprintf(1, format, d.folder, d.protocol)
    end
    fprintf(1, '\n\n');

    if nargout > 0
        out = datasets;
    end
       