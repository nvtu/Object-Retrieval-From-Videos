function tf_idf(loadDir, saveDir, fileName, saveName)

	% % loadDir = directory to load quantize file 
	% % saveDir = directory to save quantize file after using tfidf
	% % fileName = name of quantize file needed to load
	% % saveName = name of quantize file after using tfidf needed to save

    q=load(fullfile(loadDir,fileName));
    
    qtize=TFIDF(q.qtize);

    save(fullfile(saveDir,saveName),'qtize');
end

%********CALCULATE TFIDF*****************
function tfidf = TFIDF(vecBOW)
    numDocument = size(vecBOW,2);
    idf = IDF(numDocument, vecBOW);
    idf = repmat(idf,1,numDocument);
    tfidf = vecBOW.*idf;
end

%********************IDF CALCULATING************************
%        N: number of documents                            *
%        k: number of words                                *
%        len: number of images                             *
%        src: term frequency (vector BOW)                  *
%***********************************************************

function idf = IDF(numDocument,src)
    N = numDocument;
    k = size(src,1);
    idf = zeros(k,1);
    for i = 1 : k
        cnt = sum(src(i,:)~=0);
        rank = log(N/(1+cnt));
        idf(i) = rank;
    end
end
%******************************************