function asymmetric_comparison(loadDir, saveDir, fileName, queryName, saveName)

	% % loadDir = directory to load quantize file 
	% % saveDir = directory to save quantize file after using tfidf
	% % fileName = name of quantize file needed to load
	% % queryName = name of quantize file of extracted frames from video
	% % saveName = name of quantize file after using tfidf needed to save

    qtize=load(fullfile(loadDir,fileName));
    qtize=qtize.qtize;
    db_bow = double(qtize);
    ivf = BuildInvFile([],db_bow,0,false);

    qtize=load(fullfile(loadDir,queryName));
    qtize=qtize.qtize;
    qtize = double(qtize);

    dist = comp_dist(ivf,qtize,[],'autoasym_ivf_0.5',false);

    save(fullfile(saveDir,saveName),'dist');
end