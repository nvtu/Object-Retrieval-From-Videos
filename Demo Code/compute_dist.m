function dist = compute_dist(query_bow, db_bow, para)
%function dist = compute_dist(query_bow, query_fg_bow, db_bow, para)
%if strcmp(para{1},'qonly')
%    db_bow = db_bow(query_bow~=0,:);
%    query_bow = query_bow(query_bow~=0);
%    tic;
%    switch para{2}
%    case 'l1'
%        for i=1:size(db_bow,2)
%            db_bow(:,i) = db_bow(:,i)./norm(db_bow(:,i),1);
%        end
%        query_bow = query_bow./norm(query_bow,1);
%    case 'l1root'
%        for i=1:size(db_bow,2)
%            db_bow(:,i) = db_bow(:,i)./norm(db_bow(:,i),1);
%        end
%        db_bow = sqrt(db_bow);
%        query_bow = sqrt(query_bow./norm(query_bow,1));
%    case 'l2'
%        for i=1:size(db_bow,2)
%            db_bow(:,i) = db_bow(:,i)./norm(db_bow(:,i),2);
%        end
%        query_bow = query_bow./norm(query_bow,2);
%    otherwise
%    end
%    fprintf('normalization time %.0f\n',toc);
%    dist_method = para{3};
%else
%    dist_method = para{1};
%end
dist_method = para{1};
switch dist_method 
    case {'l1','l2'}
        dist=vl_alldist2(db_bow,query_bow,dist_method);
    case 'intersect'
        dist=-sum(min(repmat(query_bow,1,size(db_bow,2)),db_bow),1);
    case 'interplus'
        inter=sum(min(repmat(query_bow,1,size(db_bow,2)),db_bow),1);
        coe = str2double(para{end});
        summation = sum(query_bow)*coe+sum(db_bow,1);
        dist = -(coe+1)*inter./summation;
    case 'intermulti'
        inter=sum(min(repmat(query_bow,1,size(db_bow,2)),db_bow),1);
        multi = sum(query_bow)*sum(db_bow,1);
        dist = -(inter.*inter)./multi;
    case 'asymetric'
        dist = zeros(size(db_bow,2),1);
        coe1 = str2double(para{end-1});
        coe2 = str2double(para{end});
        for i=1:size(db_bow,2)
            big_db_idx = db_bow(:,i) > query_bow;
            big_db_diff = db_bow(big_db_idx,i)-query_bow(big_db_idx);
            big_fg_idx = query_fg_bow>db_bow(:,i);
            big_fg_diff = query_fg_bow(big_fg_idx)-db_bow(big_fg_idx,i);
            dist(i) = sum(big_db_diff)*coe1+sum(big_fg_diff)*coe2;
        end
    case 'autoasym'
        sum_qbig = zeros(1,size(db_bow,2));
        sum_qsma = zeros(1,size(db_bow,2));
        coe1 = str2double(para{end-1});
        coe2 = str2double(para{end});
        for i=1:size(db_bow,2)
            diff1 = query_fg_bow-db_bow(:,i);
            diff2 = db_bow(:,i) - query_bow;
            sum_qbig(i) = sum(diff1(diff1>0));
            sum_qsma(i) = sum(diff2(diff2>0));
        end
        mean_qbig = mean(sum_qbig);
        mean_qsma = mean(sum_qsma);
        coe2 = coe2*mean_qsma/mean_qbig;
        dist = zeros(size(db_bow,2),1);
        for i=1:size(db_bow,2)
            big_db_idx = db_bow(:,i) > query_bow;
            big_db_diff = db_bow(big_db_idx,i)-query_bow(big_db_idx);
            big_fg_idx = query_fg_bow>db_bow(:,i);
            big_fg_diff = query_fg_bow(big_fg_idx)-db_bow(big_fg_idx,i);
            dist(i) = sum(big_db_diff)*coe1+sum(big_fg_diff)*coe2;
        end
        fprintf('mean_qbig: %.02f mean_qsma: %.02f coe1: %.02f coe2: %.02f\n',mean_qbig, mean_qsma, coe1, coe2);
end
