function norm2_data(loadDir, saveDir, fileName, saveName)
	
	% % loadDir = directory to load quantize file 
	% % saveDir = directory to save quantize file after using norm2
	% % fileName = name of quantize file needed to load
	% % saveName = name of quantize file after using norm2 needed to save

    q=load(fullfile(loadDir,fileName));
    qtize=q.qtize;
    tmp = repmat(sum(qtize.^2).^0.5,size(qtize,1),1);
    qtize = qtize./tmp;

    save(fullfile(saveDir,saveName),'qtize');
end

