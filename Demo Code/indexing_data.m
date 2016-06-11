function [ind, ind2]=indexing_data(photoDir)
    
    % % photoDir = directory of database images

    list = dir(photoDir);
    n = size(list,1);
    ind = [];
    ind2 = [];
    d = 1;
    for i = 3 : n
        ind{d} = fullfile(photoDir,list(i).name);
        ind2{d} = list(i).name;
        d = d + 1;
    end
end