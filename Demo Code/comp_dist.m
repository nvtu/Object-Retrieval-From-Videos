function distance = comp_dist(ivf, query_bow, db_bow, opt_sentence, verbose, bridges, word_ini, query_b, th)
if(nargin>6)
    sim = SearchInvFile(ivf,query_bow,1,verbose, word_ini, query_b, th);
    distance = 2-2*sim;

elseif strcmp(opt_sentence,'l1_ivf')
    %tic;
    sim = SearchInvFile(ivf,query_bow,1,verbose);
    %ivf_time = toc;
    %fprintf('Search time: ivf %.04f\n',ivf_time);
    %tic
    %dist=vl_alldist2(topic_bow{qid}{i},db_bow,'l1');
    %bf_time = toc;
    %fprintf('Search time: ivf %.04f, bf: %.04f\n',ivf_time,bf_time);
    %assert(isempty(find(abs(sim*2+dist'-2)>0.0000000000001,1,'first')));
    distance = 2-2*sim;
elseif strcmp(opt_sentence,'l2_ivf')
    %tic;
    sim = SearchInvFile(ivf,query_bow,2,verbose);
    %ivf_time = toc;
    %fprintf('Search time: ivf %.04f\n',ivf_time);
    %tic
    %dist=vl_alldist2(topic_bow{qid}{i},db_bow,'l2');
    %bf_time = toc;
    %fprintf('Search time: ivf %.04f, bf: %.04f\n',ivf_time,bf_time);
    %assert(isempty(find(abs(sim*2+dist'-2)>0.0000000000001,1,'first'))
    %);
    distance = 2-2*sim;
elseif strfind(opt_sentence,'asym_ivf')==1
    %tic;
    para = textscan(opt_sentence,'%s','delimiter', '_');
    para = str2double(para{1}(3:end));
    para = [3;para];
    distance = SearchInvFile(ivf,query_bow,para,verbose);
    %ivf_time = toc;
    %fprintf('Search time: ivf %.04f\n',ivf_time);
elseif strfind(opt_sentence,'autoasym_ivf')==1
    %tic;
    para = textscan(opt_sentence,'%s','delimiter', '_');
    para = str2double(para{1}(3:end));
    para = [4;para];
    distance = SearchInvFile(ivf,query_bow,para,verbose);
    %ivf_time = toc;
    %fprintf('Search time: ivf %.04f\n',ivf_time);
elseif strfind(opt_sentence,'l2asym_ivf')==1
    %tic;
    para = textscan(opt_sentence,'%s','delimiter', '_');
    para = str2double(para{1}(3:end));
    para = [5;para];
    distance = SearchInvFile(ivf,query_bow,para,verbose);
    %ivf_time = toc;
    %fprintf('Search time: ivf %.04f\n',ivf_time);
elseif strfind(opt_sentence,'l2autoasym_ivf')==1
    %tic;
    para = textscan(opt_sentence,'%s','delimiter', '_');
    para = str2double(para{1}(3:end));
    para = [6;para];
    distance = SearchInvFile(ivf,query_bow,para,verbose);
    %ivf_time = toc;
    %fprintf('Search time: ivf %.04f\n',ivf_time);
elseif strfind(opt_sentence,'l2autoasymnew_ivf')==1
    %tic;
    para = textscan(opt_sentence,'%s','delimiter', '_');
    para = str2double(para{1}(3:end));
    para = [7;para];
    distance = SearchInvFile(ivf,query_bow,para,verbose);
    %ivf_time = toc;
    %fprintf('Search time: ivf %.04f\n',ivf_time);
elseif strfind(opt_sentence,'l2autoasymnewsqr_ivf')==1
    %tic;
    para = textscan(opt_sentence,'%s','delimiter', '_');
    para = str2double(para{1}(3:end));
    para = [8;para];
    distance = SearchInvFile(ivf,query_bow,para,verbose);
    %ivf_time = toc;
    %fprintf('Search time: ivf %.04f\n',ivf_time);
elseif strcmp(opt_sentence,'l1_bridge_ivf')
    %tic;
    sim = SearchInvFile(ivf,query_bow,9,verbose,bridges);
    distance = 2-2*sim;
elseif strcmp(opt_sentence,'l2_bridge_ivf')
    %tic;
    sim = SearchInvFile(ivf,query_bow,10,verbose,bridges);
    distance = 2-2*sim;
else
    para = textscan(opt_sentence,'%s','delimiter', '_');
    distance = compute_dist(query_bow,db_bow,para{1});
end
