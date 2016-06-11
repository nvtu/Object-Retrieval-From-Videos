function [indq, indqn]=indexing_query(queryDir)
    
    % % queryDir = directory of extracted frames of the query videos 

    list = dir(fullfile(queryDir,'*.jpg'));
    n = size(list,1);
    d = 1;
    indq = [];
    for j = 1 : n
        indq{d} = fullfile(queryDir,list(j).name);
        indqn{d} = list(j).name;
        d = d + 1;
    end
end