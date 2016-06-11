function ivf = BuildInvFile(ivf, db_bow, db_id_offset, verbose, words_ini, clip_b)

if ~exist('ivf','var') || isempty(ivf), ivf = 0; end;
if ~exist('db_id_offset','var') || isempty(db_id_offset), db_id_offset = 0; end;
if ~exist('verbose','var') || isempty(verbose), verbose = false; end;
if ~issparse(db_bow)
    db_bow = sparse(db_bow);
end
if nargin>4
    ivf = mxBuildInvFile(ivf,db_bow,db_id_offset,verbose, int32(words_ini), clip_b);
else
    ivf = mxBuildInvFile(ivf,db_bow,db_id_offset,verbose);
end