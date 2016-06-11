function sim = SearchInvFile(ivf, query_bow, para, verbose,para2, query_b, th)
if ~exist('verbose','var') || isempty(verbose), verbose = false; end;
if ~exist('para','var') || isempty(para), para = 1; end;
if ~issparse(query_bow)
    query_bow = sparse(query_bow);
end
if nargin>5
    words_ini = para2;
    sim=mxSearchInvFile(ivf,query_bow,para,verbose,int32(words_ini), logical(query_b), int32(th));
elseif nargin==5  
    bridges = para2;  
    sim=mxSearchInvFile(ivf,query_bow,para,verbose,bridges);
elseif nargin==4
    sim=mxSearchInvFile(ivf,query_bow,para,verbose);
end

