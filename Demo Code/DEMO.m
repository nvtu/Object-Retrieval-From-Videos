% % Initialization

photoDir = 'E:\APCS Program\SC203-Scientific Research\DEMO SC203\train data (photoshop)';
queryDir = 'E:\APCS Program\SC203-Scientific Research\DEMO SC203\Query';

% % TF-IDF Initialization
tfidf_loadDir = 'E:\APCS Program\SC203-Scientific Research\DEMO SC203\quantize';
tfidf_saveDir = 'E:\APCS Program\SC203-Scientific Research\DEMO SC203\tfidf';
quantizeFileName = 'qtize_features with Hessian cluster=50k';
quantizeQueryFileName = 'qtize_features Q with Hessian cluster=50k';
tfidf_saveName = 'tfidf Hessian cluster=50k';
tfidf_querySaveName = 'tfidf Q Hessian cluster=50k';

% % Normalization Initialization
norm_loadDir = 'E:\APCS Program\SC203-Scientific Research\DEMO SC203\tfidf';
norm_saveDir = 'E:\APCS Program\SC203-Scientific Research\DEMO SC203\norm';
tfidfFileName ='tfidf Hessian cluster=50k';
tfidfQueryFileName = 'tfidf Q Hessian cluster=50k';
norm_saveName = 'norm Hessian cluster=50k';
norm_querySaveName = 'norm Q Hessian cluster=50k';

% % Asymmetric Distance Calculation Initialization
asym_loadDir = 'E:\APCS Program\SC203-Scientific Research\DEMO SC203\norm';
asym_saveDir = 'E:\APCS Program\SC203-Scientific Research\DEMO SC203\asym';
normFileName =  'norm Hessian cluster=50k';
normQueryFileName = 'norm Q Hessian cluster=50k';
asym_saveName = 'asym answer Hessian cluster=50k';


% % Answer Data Initialization
ans_loadDir = 'E:\APCS Program\SC203-Scientific Research\DEMO SC203\asym';
asymFileName = 'asym answer Hessian cluster=50k';
srcImgDir = 'E:\APCS Program\SC203-Scientific Research\DEMO SC203\train data (photoshop)';
[ind, ind2] = indexing_data(photoDir);
[indq, indqn] = indexing_query(queryDir);
fprintf('Finish index process\n');

fprintf('Finish Initialization\n');

% % MAIN PROGRAM 

% tf_idf(tfidf_loadDir, tfidf_saveDir, quantizeFileName, tfidf_saveName);
% tf_idf(tfidf_loadDir, tfidf_saveDir, quantizeQueryFileName, tfidf_querySaveName);
% fprintf('Finish TF-IDF step\n');
%  
% norm2_data(norm_loadDir, norm_saveDir, tfidfFileName, norm_saveName);
% norm2_data(norm_loadDir, norm_saveDir, tfidfQueryFileName, norm_querySaveName);
% fprintf('Finish Normalization step\n');
%  
% asymmetric_comparison(asym_loadDir, asym_saveDir, normFileName, normQueryFileName, asym_saveName);
% fprintf('Finish Asymmetric Comparison step\n');

% % Copy the image needed to query to the following path: E:\APCS Program\SC203-Scientific Research\DEMO SC203\queryImg
threshold = 10;
list_ob = dir('E:\APCS Program\SC203-Scientific Research\DEMO SC203\queryImg');
num_ob = size(list_ob,1);
for i = 3 : num_ob
	objectFileName = list_ob(i).name;
	answer_data(ans_loadDir, photoDir, asymFileName, objectFileName, threshold, ind, indq, indqn);
end