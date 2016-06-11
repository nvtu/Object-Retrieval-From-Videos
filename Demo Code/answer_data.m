function answer_data(loadDir, srcImgDir, fileName, objectFileName, threshold, ind, indq, indqn)
    
    dist=load(fullfile(loadDir,fileName));
    dist=dist.dist;
    dist = dist';
    [dist2, I] = sort(dist,1);
   
    for i = 1:size(dist,2)
        if (strcmp(ind{i},fullfile(srcImgDir,objectFileName)))
            figure;
            imshow(ind{i});

            % % Convert file names to real time second        
            name=indqn(I(:,i))';
            nameList=name(:,1);
            second=[];
            for j = 1 : 2837
                name=nameList{j};
                idx=str2double(name(2:(length(name)-4)));
                if (idx>2708)
                    second=[second,int64(idx-2709+1200)];
                else
                    second=[second,int64((idx*7+1)/25)];
                end
            end

            for j = 1:threshold
                figure('Name',strcat(num2str(second(j)/3600),':',num2str(mod(second(j),3600)/60),':',num2str(mod(mod(second(j),3600),60))),'NumberTitle','off');
                imshow(indq{I(j,i)}); 
            end
            input('press any key to continue');
            close all;

            break;
        end
    end
end