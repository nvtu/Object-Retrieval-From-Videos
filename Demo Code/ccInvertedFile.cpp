#include <cstring>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <iterator>
#include <vector>
#include <map>
#include <assert.h>

#include "ccInvertedFile.hpp"

//fill the inverted file with input counts
//
// data     - the input data, with one data vector per input consisting of
//            all the word labels for its tokens
void ivBuildInvFile(ivFile* ivf, const mxArray* db, size_t docOffset, bool verbose)
{
    ivf->nwords = mxGetM(db);
    ivf->ndocs = mxGetN(db);

    //allocate vectors
    ivf->words.resize(ivf->nwords);
    mxClassID classID = mxGetClassID(db);
    assert(classID == mxDOUBLE_CLASS);
    assert(mxIsSparse(db));
    double* db_pr = (double*) mxGetData(db);
    size_t* db_jc = (size_t*) mxGetJc(db);
    size_t* db_ir = (size_t*) mxGetIr(db);
    unsigned int progress = 100;
    unsigned int percentage = ivf->ndocs/progress;
    unsigned int verb_thre = 0;
    double nz = mxGetNzmax(db);
    ivf->sparsity = nz/ivf->nwords/ivf->ndocs;
    if(verbose){
        cout<<"database sparsity: "<<ivf->sparsity<<endl;
    }
    ivf->psumofdoc = new double[ivf->ndocs];
    ivf->psqrsumofdoc = new double[ivf->ndocs];
    for(unsigned int col=0; col<ivf->ndocs; col++) {
        //get the staring index for this column
        size_t rstart = db_jc[col];
        size_t rend = db_jc[col+1];
        double& sumofdoc = ivf->psumofdoc[col];
        sumofdoc = 0;
        double& sqrsumofdoc = ivf->psqrsumofdoc[col];
        sqrsumofdoc = 0;
        if (isnormal(*(db_pr+rstart))){
            //column has something
            for (size_t r=rstart; r<rend; r++) {
                double tf = *(db_pr+r);
                ivWordDoc wordDoc;
                sumofdoc += tf;
                sqrsumofdoc += tf*tf;
                wordDoc.termFreq = tf;
                wordDoc.docID=col+docOffset;
                ivf->words[*(db_ir+r)].push_back(wordDoc);
                //cout<<"word "<<*(db_ir+r)<<": doc "<<col<<" val "<<*(db_pr+r)<<endl;
            }
        }
        if (verbose){
            if(percentage==0){
                mexPrintf("\r%d/%d,sumofdoc:%.2f",col+1,ivf->ndocs,ivf->psumofdoc[col]);
                mexEvalString("drawnow;"); // to print string immediately.
               // cout<<col+1<<"/"<<ivf->ndocs<<endl;
            }
            else if(col>=verb_thre){
                mexPrintf("\r%d%%,sumofdoc:%.2f",col/percentage,ivf->psumofdoc[col]);
                mexEvalString("drawnow;"); // to print string immediately.
                //cout<<col/percentage*progress<<"%"<<endl;
                verb_thre += percentage;
            }
        }
    }
}

// the input db here acts as a forward indexing file, can be boolean
mxArray* ivComputeSelfSim(ivFile& ivf, const mxArray* db,  bool verbose)
{
    /* Declare variable */
    mwSize m,n;
    mwSize nzmax;
    mwIndex *irs,*jcs,col,count;
    int isfull;
    double *pr,*sr;
    double tf;
    double percent_sparse;

    size_t r,rstart,rend;

    size_t nwords = ivf.nwords;
    size_t ndocs = ivf.ndocs;
    unsigned int* lut = new unsigned int[ndocs];

    unsigned int progress = 100;
    unsigned int percentage = ndocs/progress;
    unsigned int verb_thre = 0;

    if (verbose){
        mexPrintf("sparsity %.0f%%\n",ivf.sparsity*100);
        mexEvalString("drawnow;"); // to print string immediately.
    }
    percent_sparse = max(0.01,ivf.sparsity);
    nzmax=(mwSize)ceil(ndocs*ndocs*percent_sparse);
    mxArray* matSelfSim = mxCreateSparse(ndocs, ndocs, nzmax, mxREAL);
    sr  = mxGetPr(matSelfSim);
    irs = mxGetIr(matSelfSim);
    jcs = mxGetJc(matSelfSim);

    assert(mxIsSparse(db));
    double* db_pr = (double*) mxGetData(db);
    size_t* db_jc = (size_t*) mxGetJc(db);
    size_t* db_ir = (size_t*) mxGetIr(db);
    count = 0;
    for(col=0; col<ndocs; col++) {
        //get the staring index for this column
        rstart = db_jc[col];
        rend = db_jc[col+1];
        memset(lut, 0, sizeof(unsigned int)*ndocs);
        vector<double> vecSim;
        // compute the similarity between the col-th doc with all others
        if (isnormal(*(db_pr+rstart))){
            //column has something, i.e., the col-th doc is not empty
            for (r=rstart; r<rend; r++) {
                // get nonzero word entries of the col-th doc
                int wordID = *(db_ir+r);
                float termFreq = *(db_pr+r);
                for (ivWordDocIt it=ivf.words[wordID].begin(); it!=ivf.words[wordID].end(); it++)
                    // acumulate similarities for all docs in the inverted list cooresponding the word
                {
                    tf = min(termFreq,it->termFreq);
                    if (lut[it->docID]==0) // is new
                    {
                        vecSim.push_back(tf);
                        lut[it->docID]=vecSim.size();
                    }
                    else
                        vecSim[lut[it->docID]-1] += tf;
                }
            }
        }
        if (false && verbose && col==0){
#define N 4
            double big_sim[N]={0,0,0,0};
            int big_sim_idx[N];
            vector<int> zero_sim_idx;
            double itself_sim;

            mexPrintf("****** image idx %d ******\n",col);
            mexPrintf("self sim:[%d %g]\n",col, vecSim[lut[col]-1]);
            mexPrintf("zero sim idx:\n");
            int cnt = 0;
            for (int i=0; i<ndocs; i++){
                if (lut[i]==0){
                    mexPrintf("%d ",i);
                    if(++cnt%10 == 0)
                        mexPrintf("\n");
                }
                else{
                    double sim = vecSim[lut[i]-1];
                    for (int j=0;j<N;j++){
                        if (sim>big_sim[j]){
                            for (int k=N-1;k>j;k--){
                                big_sim[k]=big_sim[k-1];
                                big_sim_idx[k]=big_sim_idx[k-1];
                            }
                            big_sim[j]=sim;
                            big_sim_idx[j]=i;
                            break;
                        }
                    }
                }
            }
            mexPrintf("\ntop %d [idx,sim]:\n", N);
            cnt = 0;
            for (int i=0; i<N; i++){
                mexPrintf("[%d %g] ",big_sim_idx[i],big_sim[i]);
                if(++cnt%10 == 0)
                    mexPrintf("\n");
            }
            mexPrintf("\n");
            mexEvalString("drawnow;"); // to print string immediately.
        }

        // assign similarities of the col-th doc to the corresponding column in the sparse matrix
        jcs[col] = count;
        for (unsigned int i=0; i<ndocs; i++){
            if (lut[i] == 0)
                continue;

            // if similarity matrix is denser than the allocated sparse mat
            // then reallocate memory
            if (count>=nzmax){
                mwSize oldnzmax = nzmax;
                percent_sparse += 0.01;
                nzmax=(mwSize)ceil(ndocs*ndocs*percent_sparse);

                /* make sure nzmax increases atleast by 1 */
                if (oldnzmax == nzmax)
                    nzmax++;

                mxSetNzmax(matSelfSim, nzmax);
                mxSetPr(matSelfSim, (double*)mxRealloc(sr, nzmax*sizeof(double)));
                mxSetIr(matSelfSim, (mwIndex*)mxRealloc(irs, nzmax*sizeof(mwIndex)));

                sr  = mxGetPr(matSelfSim);
                irs = mxGetIr(matSelfSim);
            }

            // assignment
            sr[count] = vecSim[lut[i]-1];
            irs[count] = i;

            if (false && verbose && col<10 && irs[count]<10){
                mexPrintf("(%d, %d)\t%g\n",col+1,irs[count]+1,sr[count]);
            }
            count++;
        }

        if (verbose && col>=verb_thre){
            mexPrintf("\r%d%%",col/percentage);
            mexEvalString("drawnow;"); // to print string immediately.
            //cout<<col/percentage*progress<<"%"<<endl;
            verb_thre += percentage;
        }
    }
    if (verbose)
        mexEvalString("drawnow;"); // to print string immediately.
    jcs[ndocs] = count;

    delete []lut;
    return matSelfSim;
}

mxArray* ivComputeSelfSim_l1(ivFile& ivf, const mxArray* db,  bool verbose)
{
    ivComputeSelfSim(ivf, db, verbose);
}

// the input db here acts as a forward indexing file, can be boolean
mxArray* ivComputeSelfSim_l2(ivFile& ivf, const mxArray* db,  bool verbose)
{
    /* Declare variable */
    mwSize m,n;
    mwSize nzmax;
    mwIndex *irs,*jcs,col,count;
    int isfull;
    double *pr,*sr;
    double tf;
    double percent_sparse;

    size_t r,rstart,rend;

    size_t nwords = ivf.nwords;
    size_t ndocs = ivf.ndocs;
    unsigned int* lut = new unsigned int[ndocs];

    unsigned int progress = 100;
    unsigned int percentage = ndocs/progress;
    unsigned int verb_thre = 0;

    if (verbose){
        mexPrintf("sparsity %.0f%%\n",ivf.sparsity*100);
        mexEvalString("drawnow;"); // to print string immediately.
    }
    percent_sparse = max(0.01,ivf.sparsity);
    nzmax=(mwSize)ceil(ndocs*ndocs*percent_sparse);
    mxArray* matSelfSim = mxCreateSparse(ndocs, ndocs, nzmax, mxREAL);
    sr  = mxGetPr(matSelfSim);
    irs = mxGetIr(matSelfSim);
    jcs = mxGetJc(matSelfSim);

    assert(mxIsSparse(db));
    double* db_pr = (double*) mxGetData(db);
    size_t* db_jc = (size_t*) mxGetJc(db);
    size_t* db_ir = (size_t*) mxGetIr(db);
    count = 0;
    for(col=0; col<ndocs; col++) {
        //get the staring index for this column
        rstart = db_jc[col];
        rend = db_jc[col+1];
        memset(lut, 0, sizeof(unsigned int)*ndocs);
        vector<double> vecSim;
        // compute the similarity between the col-th doc with all others
        if (isnormal(*(db_pr+rstart))){
            //column has something, i.e., the col-th doc is not empty
            for (r=rstart; r<rend; r++) {
                // get nonzero word entries of the col-th doc
                int wordID = *(db_ir+r);
                float termFreq = *(db_pr+r);
                for (ivWordDocIt it=ivf.words[wordID].begin(); it!=ivf.words[wordID].end(); it++)
                    // acumulate similarities for all docs in the inverted list cooresponding the word
                {
                    tf = termFreq * it->termFreq;
                    if (lut[it->docID]==0) // is new
                    {
                        vecSim.push_back(tf);
                        lut[it->docID]=vecSim.size();
                    }
                    else
                        vecSim[lut[it->docID]-1] += tf;
                }
            }
        }
        if (false && verbose && col==0){
#define N 4
            double big_sim[N]={0,0,0,0};
            int big_sim_idx[N];
            vector<int> zero_sim_idx;
            double itself_sim;

            mexPrintf("****** image idx %d ******\n",col);
            mexPrintf("self sim:[%d %g]\n",col, vecSim[lut[col]-1]);
            mexPrintf("zero sim idx:\n");
            int cnt = 0;
            for (int i=0; i<ndocs; i++){
                if (lut[i]==0){
                    mexPrintf("%d ",i);
                    if(++cnt%10 == 0)
                        mexPrintf("\n");
                }
                else{
                    double sim = vecSim[lut[i]-1];
                    for (int j=0;j<N;j++){
                        if (sim>big_sim[j]){
                            for (int k=N-1;k>j;k--){
                                big_sim[k]=big_sim[k-1];
                                big_sim_idx[k]=big_sim_idx[k-1];
                            }
                            big_sim[j]=sim;
                            big_sim_idx[j]=i;
                            break;
                        }
                    }
                }
            }
            mexPrintf("\ntop %d [idx,sim]:\n", N);
            cnt = 0;
            for (int i=0; i<N; i++){
                mexPrintf("[%d %g] ",big_sim_idx[i],big_sim[i]);
                if(++cnt%10 == 0)
                    mexPrintf("\n");
            }
            mexPrintf("\n");
            mexEvalString("drawnow;"); // to print string immediately.
        }

        // assign similarities of the col-th doc to the corresponding column in the sparse matrix
        jcs[col] = count;
        for (unsigned int i=0; i<ndocs; i++){
            if (lut[i] == 0)
                continue;

            // if similarity matrix is denser than the allocated sparse mat
            // then reallocate memory
            if (count>=nzmax){
                mwSize oldnzmax = nzmax;
                percent_sparse += 0.01;
                nzmax=(mwSize)ceil(ndocs*ndocs*percent_sparse);

                /* make sure nzmax increases atleast by 1 */
                if (oldnzmax == nzmax)
                    nzmax++;

                mxSetNzmax(matSelfSim, nzmax);
                mxSetPr(matSelfSim, (double*)mxRealloc(sr, nzmax*sizeof(double)));
                mxSetIr(matSelfSim, (mwIndex*)mxRealloc(irs, nzmax*sizeof(mwIndex)));

                sr  = mxGetPr(matSelfSim);
                irs = mxGetIr(matSelfSim);
            }

            // assignment
            sr[count] = vecSim[lut[i]-1];
            irs[count] = i;

            if (false && verbose && col<10 && irs[count]<10){
                mexPrintf("(%d, %d)\t%g\n",col+1,irs[count]+1,sr[count]);
            }
            count++;
        }

        if (verbose && col>=verb_thre){
            mexPrintf("\r%d%%",col/percentage);
            mexEvalString("drawnow;"); // to print string immediately.
            //cout<<col/percentage*progress<<"%"<<endl;
            verb_thre += percentage;
        }
    }
    if (verbose)
        mexEvalString("drawnow;"); // to print string immediately.
    jcs[ndocs] = count;

    delete []lut;
    return matSelfSim;
}

mxArray* ivSearchInvFile_l1(ivFile& ivf, const mxArray* queries, bool verbose)
{
    size_t nwords, nqueries;
    nwords = mxGetM(queries);
    nqueries = mxGetN(queries);
    assert(nwords == ivf.nwords);
    mxArray* sim = mxCreateNumericMatrix(ivf.ndocs, nqueries, mxSINGLE_CLASS, mxREAL);
    float *psim = (float*)mxGetPr(sim);
    mxClassID classID = mxGetClassID(queries);
    assert(classID == mxDOUBLE_CLASS);
    assert(mxIsSparse(queries));
    unsigned int progress = 100;
    unsigned int percentage = nqueries/progress;
    unsigned int verb_thre = 0;
    if(verbose){
        double nz = mxGetNzmax(queries);
        cout<<"query sparsity: "<<nz/nwords/nqueries<<endl;
    }
    double* query_pr = (double*) mxGetData(queries);
    size_t* query_jc = (size_t*) mxGetJc(queries);
    size_t* query_ir = (size_t*) mxGetIr(queries);
    for(size_t col=0; col<nqueries; col++)
    {
        //get the staring index for this column
        size_t rstart = query_jc[col];
        size_t rend = query_jc[col+1];
        if (isnormal(*(query_pr+rstart))){
            //column has something
            for (size_t r=rstart; r<rend; r++)
            {
                float tf = float(*(query_pr+r));
                size_t word_id = *(query_ir+r);
                for (ivWordDocIt it=ivf.words[word_id].begin(); it!=ivf.words[word_id].end(); it++)
                    psim[it->docID] += min(it->termFreq,tf);
            }
        }
        psim += ivf.ndocs;
        if (verbose && nqueries>1){
            if(percentage==0){
                mexPrintf("\r%d/%d",col+1,nqueries);
                mexEvalString("drawnow;"); // to print string immediately.
                //cout<<col+1<<"/"<<nqueries<<endl;
            }
            else if(col>=verb_thre){
                mexPrintf("\r%d%%",col/percentage);
                mexEvalString("drawnow;"); // to print string immediately.
                //cout<<col/percentage<<"%"<<endl;
                verb_thre += percentage;
            }
        }
    }
    return sim;
}

mxArray* ivSearchInvFile_bridge_l1(ivFile& ivf, const mxArray* queries, const mxArray* query_bridges, bool verbose)
{
    size_t nwords, nqueries;
    nwords = mxGetM(queries);
    nqueries = mxGetN(queries);
    assert(nwords == ivf.nwords);
    mxArray* sim = mxCreateNumericMatrix(ivf.ndocs, nqueries, mxSINGLE_CLASS, mxREAL);
    float *psim = (float*)mxGetPr(sim);
    mxClassID classID = mxGetClassID(queries);
    assert(classID == mxDOUBLE_CLASS);
    assert(mxIsSparse(queries));
    unsigned int progress = 100;
    unsigned int percentage = nqueries/progress;
    unsigned int verb_thre = 0;
    if(verbose){
        double nz = mxGetNzmax(queries);
        cout<<"query sparsity: "<<nz/nwords/nqueries<<endl;
    }
    double* query_pr = (double*) mxGetData(queries);
    size_t* query_jc = (size_t*) mxGetJc(queries);
    size_t* query_ir = (size_t*) mxGetIr(queries);

    if (!mxIsCell(query_bridges))
        mexErrMsgTxt("Need a cell array");
    for(size_t col=0; col<nqueries; col++)
    {
        //get the query bridges
        mxArray *query_bridge = mxGetCell(query_bridges,col);
        if (!mxIsCell(query_bridge))
            mexErrMsgTxt("Need a cell array");

        int numel = mxGetNumberOfElements(query_bridge);
        vector<vector<int> > vv_bridge;
        if(verbose) cout<<numel<<" bridges:"<<endl;
        for (int i=0; i<numel; i++)
        {
            mxArray *pArrayIn = mxGetCell(query_bridge,i);
            double *pDataIn = (double *)mxGetData(pArrayIn);
            int num_words =  mxGetM(pArrayIn)*mxGetN(pArrayIn); //Number of rows
            vector<int> vec_bridge;
            for (int j=0; j<num_words; j++){
                vec_bridge.push_back(int(pDataIn[j]-1)); //-1 because of matlab start from 1
                if(verbose) cout<<vec_bridge[j]+1<<" ";
            }
            if(verbose) cout<<endl;
            vv_bridge.push_back(vec_bridge);
        }


        //get the staring index for this column
        size_t rstart = query_jc[col];
        size_t rend = query_jc[col+1];
        vector<bool> vec_comp_bridge;
        vec_comp_bridge.assign(numel,false);
        if (isnormal(*(query_pr+rstart))){
//            if(verbose){
//                for (size_t r=rstart; r<rend; r++)
//                    cout<<*(query_ir+r)<<" ";
//                cout<<endl;
//            }

            //column has something
            for (size_t r=rstart; r<rend; r++)
            {
                size_t word_id = *(query_ir+r);
                bool bridged = false;
                int i=0;
                for(i=0; i<numel; i++){
                    for(int j=0; j<vv_bridge[i].size(); j++){
                        if (word_id == vv_bridge[i][j]){
                            bridged = true;
                            break;
                        }
                    }
                    if(bridged)
                        break;
                }

                if(!bridged){
                    float tf = float(*(query_pr+r));
                    for (ivWordDocIt it=ivf.words[word_id].begin(); it!=ivf.words[word_id].end(); it++)
                        psim[it->docID] += min(it->termFreq,tf);
                }
                else{
                    // compute each bridge only once.
                    if (vec_comp_bridge[i])
                        continue;

                    float tf = 0;
                    int dup=0;
                    map<unsigned int, float> map_id_tf;
                    for(int j=0; j<vv_bridge[i].size(); j++){
                        word_id = vv_bridge[i][j];
                        // compute accumulated tf on this bridge for the query.
                        for (size_t q=rstart; q<rend; q++){
                            if(word_id == *(query_ir+q)){
                                tf += float(*(query_pr+q));
                                break;
                            }
                        }
                        // compute accumulated tf on this bridge for all database items.
                        for (ivWordDocIt it=ivf.words[word_id].begin(); it!=ivf.words[word_id].end(); it++){
                            if (map_id_tf.find(it->docID) == map_id_tf.end())
                                map_id_tf.insert(pair<unsigned int,float>(it->docID,it->termFreq));
                            else{
                                map_id_tf[it->docID] = map_id_tf[it->docID] + it->termFreq;
                                dup++;
                            }
                        }
                    }

                    for (map<unsigned int,float>::iterator mit=map_id_tf.begin(); mit!=map_id_tf.end(); mit++)
                        psim[mit->first] += min(mit->second,tf);

                    if (verbose) cout<<"map size "<<map_id_tf.size()<<"; dup:"<<dup<<endl;
                    vec_comp_bridge[i] = true;
                }
            }
        }
        if (verbose){
            for (size_t n=0;n<vec_comp_bridge.size();n++)
                cout<<vec_comp_bridge[n]<<" ";
            cout<<endl<<"-----------"<<endl;
        }
        psim += ivf.ndocs;
        if (verbose && nqueries>1){
            if(percentage==0){
                mexPrintf("\r%d/%d",col+1,nqueries);
                mexEvalString("drawnow;"); // to print string immediately.
                //cout<<col+1<<"/"<<nqueries<<endl;
            }
            else if(col>=verb_thre){
                mexPrintf("\r%d%%",col/percentage);
                mexEvalString("drawnow;"); // to print string immediately.
                //cout<<col/percentage<<"%"<<endl;
                verb_thre += percentage;
            }
        }
    }
    return sim;
}

mxArray* ivSearchInvFile_l2(ivFile& ivf, const mxArray* queries, bool verbose)
{
    size_t nwords, nqueries;
    nwords = mxGetM(queries);
    nqueries = mxGetN(queries);
    assert(nwords == ivf.nwords);
    mxArray* sim = mxCreateNumericMatrix(ivf.ndocs, nqueries, mxSINGLE_CLASS, mxREAL);
    float *psim = (float*)mxGetPr(sim);
    mxClassID classID = mxGetClassID(queries);
    assert(classID == mxDOUBLE_CLASS);
    assert(mxIsSparse(queries));
    unsigned int progress = 100;
    unsigned int percentage = nqueries/progress;
    unsigned int verb_thre = 0;
    if(verbose){
        double nz = mxGetNzmax(queries);
        cout<<"query sparsity: "<<nz/nwords/nqueries<<endl;
    }
    double* query_pr = (double*) mxGetData(queries);
    size_t* query_jc = (size_t*) mxGetJc(queries);
    size_t* query_ir = (size_t*) mxGetIr(queries);
    for(size_t col=0; col<nqueries; col++)
    {
        //get the staring index for this column
        size_t rstart = query_jc[col];
        size_t rend = query_jc[col+1];
        if (isnormal(*(query_pr+rstart))){
            //column has something
            for (size_t r=rstart; r<rend; r++)
            {
                float tf = float(*(query_pr+r));
                size_t word_id = *(query_ir+r);
                for (ivWordDocIt it=ivf.words[word_id].begin(); it!=ivf.words[word_id].end(); it++)
                    psim[it->docID] += it->termFreq * tf;
            }
        }
        psim += ivf.ndocs;
        if (verbose && nqueries>1){
            if(percentage==0){
                mexPrintf("\r%d/%d",col+1,nqueries);
                mexEvalString("drawnow;"); // to print string immediately.
                //cout<<col+1<<"/"<<nqueries<<endl;
            }
            else if(col>=verb_thre){
                mexPrintf("\r%d%%",col/percentage);
                mexEvalString("drawnow;"); // to print string immediately.
                //cout<<col/percentage<<"%"<<endl;
                verb_thre += percentage;
            }
        }
    }
    return sim;
}

mxArray* ivSearchInvFile_bridge_l2(ivFile& ivf, const mxArray* queries, const mxArray* query_bridges, bool verbose)
{
    size_t nwords, nqueries;
    nwords = mxGetM(queries);
    nqueries = mxGetN(queries);
    assert(nwords == ivf.nwords);
    mxArray* sim = mxCreateNumericMatrix(ivf.ndocs, nqueries, mxSINGLE_CLASS, mxREAL);
    float *psim = (float*)mxGetPr(sim);
    mxClassID classID = mxGetClassID(queries);
    assert(classID == mxDOUBLE_CLASS);
    assert(mxIsSparse(queries));
    unsigned int progress = 100;
    unsigned int percentage = nqueries/progress;
    unsigned int verb_thre = 0;
    if(verbose){
        double nz = mxGetNzmax(queries);
        cout<<"query sparsity: "<<nz/nwords/nqueries<<endl;
    }
    double* query_pr = (double*) mxGetData(queries);
    size_t* query_jc = (size_t*) mxGetJc(queries);
    size_t* query_ir = (size_t*) mxGetIr(queries);

    if (!mxIsCell(query_bridges))
        mexErrMsgTxt("Need a cell array");
    for(size_t col=0; col<nqueries; col++)
    {
        //get the query bridges
        mxArray *query_bridge = mxGetCell(query_bridges,col);
        if (!mxIsCell(query_bridge))
            mexErrMsgTxt("Need a cell array");

        int numel = mxGetNumberOfElements(query_bridge);
        vector<vector<int> > vv_bridge;
        if(verbose) cout<<numel<<" bridges:"<<endl;
        for (int i=0; i<numel; i++)
        {
            mxArray *pArrayIn = mxGetCell(query_bridge,i);
            double *pDataIn = (double *)mxGetData(pArrayIn);
            int num_words =  mxGetM(pArrayIn)*mxGetN(pArrayIn); //Number of rows
            vector<int> vec_bridge;
            for (int j=0; j<num_words; j++){
                vec_bridge.push_back(int(pDataIn[j]-1)); //-1 because of matlab start from 1
                if(verbose) cout<<vec_bridge[j]+1<<" ";
            }
            if(verbose) cout<<endl;
            vv_bridge.push_back(vec_bridge);
        }


        //get the staring index for this column
        size_t rstart = query_jc[col];
        size_t rend = query_jc[col+1];
        vector<bool> vec_comp_bridge;
        vec_comp_bridge.assign(numel,false);
        if (isnormal(*(query_pr+rstart))){
            //column has something
            for (size_t r=rstart; r<rend; r++)
            {
                size_t word_id = *(query_ir+r);
                bool bridged = false;
                int i=0;
                for(i=0; i<numel; i++){
                    for(int j=0; j<vv_bridge[i].size(); j++){
                        if (word_id == vv_bridge[i][j]){
                            bridged = true;
                            break;
                        }
                    }
                    if(bridged)
                        break;
                }

                if(!bridged){
                    float tf = float(*(query_pr+r));
                    for (ivWordDocIt it=ivf.words[word_id].begin(); it!=ivf.words[word_id].end(); it++)
                        psim[it->docID] += it->termFreq * tf;
                }
                else{
                    // compute each bridge only once.
                    if (vec_comp_bridge[i])
                        continue;

                    float tf = 0;
                    int dup=0;
                    map<unsigned int, float> map_id_tf;
                    for(int j=0; j<vv_bridge[i].size(); j++){
                        word_id = vv_bridge[i][j];
                        // compute accumulated tf on this bridge for the query.
                        for (size_t q=rstart; q<rend; q++){
                            if(word_id == *(query_ir+q)){
                                float tf_temp = float(*(query_pr+q));
                                tf += tf_temp*tf_temp;
                                break;
                            }
                        }
                        // compute accumulated tf on this bridge for all database items.
                        for (ivWordDocIt it=ivf.words[word_id].begin(); it!=ivf.words[word_id].end(); it++){
                            if (map_id_tf.find(it->docID) == map_id_tf.end())
                                map_id_tf.insert(pair<unsigned int,float>(it->docID,it->termFreq*it->termFreq));
                            else{
                                map_id_tf[it->docID] = map_id_tf[it->docID] + it->termFreq*it->termFreq;
                                dup++;
                            }
                        }
                    }

                    for (map<unsigned int,float>::iterator mit=map_id_tf.begin(); mit!=map_id_tf.end(); mit++)
                        psim[mit->first] += sqrt(mit->second)*sqrt(tf);

                    if (verbose) cout<<"map size "<<map_id_tf.size()<<"; dup:"<<dup<<endl;
                    vec_comp_bridge[i] = true;
                }
            }
        }
        if (verbose){
            for (size_t n=0;n<vec_comp_bridge.size();n++)
                cout<<vec_comp_bridge[n]<<" ";
            cout<<endl<<"-----------"<<endl;
        }
        psim += ivf.ndocs;
        if (verbose && nqueries>1){
            if(percentage==0){
                mexPrintf("\r%d/%d",col+1,nqueries);
                mexEvalString("drawnow;"); // to print string immediately.
                //cout<<col+1<<"/"<<nqueries<<endl;
            }
            else if(col>=verb_thre){
                mexPrintf("\r%d%%",col/percentage);
                mexEvalString("drawnow;"); // to print string immediately.
                //cout<<col/percentage<<"%"<<endl;
                verb_thre += percentage;
            }
        }
    }
    return sim;
}

mxArray* ivSearchInvFile_asym(ivFile& ivf, const mxArray* queries, double w1, double w2, bool verbose)
{
    size_t nwords, nqueries;
    nwords = mxGetM(queries);
    nqueries = mxGetN(queries);
    assert(nwords == ivf.nwords);
    mxArray* dist = mxCreateNumericMatrix(ivf.ndocs, nqueries, mxSINGLE_CLASS, mxREAL);
    float *pdist = (float*)mxGetPr(dist);
    mxClassID classID = mxGetClassID(queries);
    assert(classID == mxDOUBLE_CLASS);
    assert(mxIsSparse(queries));
    unsigned int progress = 100;
    unsigned int percentage = nqueries/progress;
    unsigned int verb_thre = 0;
    if(verbose){
        double nz = mxGetNzmax(queries);
        cout<<"query sparsity: "<<nz/nwords/nqueries<<endl;
    }
    double* query_pr = (double*) mxGetData(queries);
    size_t* query_jc = (size_t*) mxGetJc(queries);
    size_t* query_ir = (size_t*) mxGetIr(queries);
    float* psumofquery = new float[nqueries];
    memset(psumofquery, 0, sizeof(float)*nqueries);
    for(size_t col=0; col<nqueries; col++)
    {
        //get the staring index for this column
        size_t rstart = query_jc[col];
        size_t rend = query_jc[col+1];
        if (isnormal(*(query_pr+rstart))){
            //column has something
            for (size_t r=rstart; r<rend; r++)
            {
                float tf = float(*(query_pr+r));
                size_t word_id = *(query_ir+r);
                psumofquery[col] += tf;
                for (ivWordDocIt it=ivf.words[word_id].begin(); it!=ivf.words[word_id].end(); it++)
                    pdist[it->docID] += min(it->termFreq,tf);
            }
        }
        pdist += ivf.ndocs;
        if (verbose){
            if (nqueries>1){
                if(percentage==0){
                    mexPrintf("\r%d/%d, sumofquery:%.2f",col+1,nqueries,psumofquery[col]);
                    mexEvalString("drawnow;"); // to print string immediately.
                    //cout<<col+1<<"/"<<nqueries<<endl;
                }
                else if(col>=verb_thre){
                    mexPrintf("\r%d%%, sumofquery:%.2f",col/percentage,psumofquery[col]);
                    mexEvalString("drawnow;"); // to print string immediately.
                    //cout<<col/percentage<<"%"<<endl;
                    verb_thre += percentage;
                }
            }
            else
                cout<<"sum of query: "<<psumofquery[col]<<endl;
        }
    }
    pdist = (float*)mxGetPr(dist);
    for(size_t col=0; col<nqueries; col++){
        for(size_t r=0; r<ivf.ndocs; r++)
            pdist[r] = (psumofquery[col]-pdist[r])*w1 + (ivf.psumofdoc[r]-pdist[r])*w2;
        pdist += ivf.ndocs;
    }
    delete []psumofquery;
    return dist;
}

mxArray* ivSearchInvFile_l2asym(ivFile& ivf, const mxArray* queries, double w1, double w2, bool verbose)
{
    if(verbose) cout<<"l2asym"<<endl;
    size_t nwords, nqueries;
    nwords = mxGetM(queries);
    nqueries = mxGetN(queries);
    assert(nwords == ivf.nwords);
    mxClassID classID = mxGetClassID(queries);
    assert(classID == mxDOUBLE_CLASS);
    assert(mxIsSparse(queries));
    unsigned int progress = 100;
    unsigned int percentage = nqueries/progress;
    unsigned int verb_thre = 0;
    double avgsqrsumofdoc = 0;
    if(verbose){
        double nz = mxGetNzmax(queries);
        cout<<"query sparsity: "<<nz/nwords/nqueries<<endl;

        for(int i=0; i<ivf.ndocs; i++)
            avgsqrsumofdoc += ivf.psqrsumofdoc[i];

        avgsqrsumofdoc /= ivf.ndocs;
    }
    double* query_pr = (double*) mxGetData(queries);
    size_t* query_jc = (size_t*) mxGetJc(queries);
    size_t* query_ir = (size_t*) mxGetIr(queries);
    double* dist1 = new double[nqueries*ivf.ndocs];
    double* dist2 = new double[nqueries*ivf.ndocs];
    double* pdist1, *pdist2;
    pdist1 = dist1; pdist2 = dist2;
    for (size_t i=0; i<nqueries; i++){
        memcpy(pdist2,ivf.psqrsumofdoc,sizeof(double)*ivf.ndocs);

        //get the staring index for this column
        size_t rstart = query_jc[i];
        size_t rend = query_jc[i+1];
        double sqrsumofquery = 0;
        if (isnormal(*(query_pr+rstart))){
            //column has something
            for (size_t r=rstart; r<rend; r++)
            {
                float tf = float(*(query_pr+r));
                sqrsumofquery += tf*tf;
            }
        }
        for(size_t j=0; j<ivf.ndocs; j++)
            pdist1[j]=sqrsumofquery;

        if(verbose)
            cout<<"avgsqrsumofdoc/sqrsumofquery: "<<avgsqrsumofdoc/sqrsumofquery<<endl;

        pdist1 += ivf.ndocs;
        pdist2 += ivf.ndocs;
    }
    pdist1 = dist1; pdist2 = dist2;
    for(size_t col=0; col<nqueries; col++)
    {
        //get the staring index for this column
        size_t rstart = query_jc[col];
        size_t rend = query_jc[col+1];
        if (isnormal(*(query_pr+rstart))){
            //column has something
            for (size_t r=rstart; r<rend; r++)
            {
                float tf = float(*(query_pr+r));
                size_t word_id = *(query_ir+r);
                for (ivWordDocIt it=ivf.words[word_id].begin(); it!=ivf.words[word_id].end(); it++){
                    float inter = min(tf,it->termFreq);
                    float sqrinter = inter*inter;
                    pdist1[it->docID] += (-2*tf*inter + sqrinter);
                    pdist2[it->docID] += (-2*it->termFreq*inter + sqrinter);
                }
            }
        }
        pdist1 += ivf.ndocs;
        pdist2 += ivf.ndocs;
    }
    mxArray* dist = mxCreateNumericMatrix(ivf.ndocs, nqueries, mxSINGLE_CLASS, mxREAL);
    float *pdist = (float*)mxGetPr(dist);
    pdist1 = dist1; pdist2 = dist2;
    for(size_t col=0; col<nqueries; col++){
        double similarity_sum = 0;
        double complexity_sum = 0;
        for(size_t r=0; r<ivf.ndocs; r++){
            pdist1[r] = sqrt(pdist1[r]);
            pdist2[r] = sqrt(pdist2[r]);
            similarity_sum += pdist1[r];
            complexity_sum += pdist2[r];
        }
        if(verbose)
            cout<<"complexity "<<complexity_sum<<" similarity "<<similarity_sum<<" ratio "<<complexity_sum/(similarity_sum+0.0001)<<endl;
        for(size_t r=0; r<ivf.ndocs; r++)
            pdist[r] = w1*pdist1[r] + w2*pdist2[r];
        pdist += ivf.ndocs;
        pdist1 += ivf.ndocs;
        pdist2 += ivf.ndocs;
    }
    delete []dist1;
    delete []dist2;
    return dist;
}

// correct the bug of zero intersection bows being ranked high.
mxArray* ivSearchInvFile_autoasym(ivFile& ivf, const mxArray* queries, double w, bool verbose)
{
    size_t nwords, nqueries;
    nwords = mxGetM(queries);
    nqueries = mxGetN(queries);
    assert(nwords == ivf.nwords);
    mxArray* dist = mxCreateNumericMatrix(ivf.ndocs, nqueries, mxSINGLE_CLASS, mxREAL);
    float *pdist = (float*)mxGetPr(dist);
    mxClassID classID = mxGetClassID(queries);
    assert(classID == mxDOUBLE_CLASS);
    assert(mxIsSparse(queries));
    unsigned int progress = 100;
    unsigned int percentage = nqueries/progress;
    unsigned int verb_thre = 0;
    if(false){
        double nz = mxGetNzmax(queries);
        cout<<"query sparsity: "<<nz/nwords/nqueries<<endl;
    }
    double* query_pr = (double*) mxGetData(queries);
    size_t* query_jc = (size_t*) mxGetJc(queries);
    size_t* query_ir = (size_t*) mxGetIr(queries);
    float* psumofquery = new float[nqueries];
    memset(psumofquery, 0, sizeof(float)*nqueries);
    for(size_t col=0; col<nqueries; col++)
    {
        //get the staring index for this column
        size_t rstart = query_jc[col];
        size_t rend = query_jc[col+1];
        if (isnormal(*(query_pr+rstart))){
            //column has something
            for (size_t r=rstart; r<rend; r++)
            {
                float tf = float(*(query_pr+r));
                size_t word_id = *(query_ir+r);
                psumofquery[col] += tf;
                for (ivWordDocIt it=ivf.words[word_id].begin(); it!=ivf.words[word_id].end(); it++)
                    pdist[it->docID] += min(it->termFreq,tf);
            }
        }
        pdist += ivf.ndocs;
        if (false){
            if (nqueries>1){
                if(percentage==0){
                    mexPrintf("\r%d/%d, sumofquery:%.2f",col+1,nqueries,psumofquery[col]);
                    mexEvalString("drawnow;"); // to print string immediately.
                    //cout<<col+1<<"/"<<nqueries<<endl;
                }
                else if(col>=verb_thre){
                    mexPrintf("\r%d%%, sumofquery:%.2f",col/percentage,psumofquery[col]);
                    mexEvalString("drawnow;"); // to print string immediately.
                    //cout<<col/percentage<<"%"<<endl;
                    verb_thre += percentage;
                }
            }
            else
                cout<<"sum of query: "<<psumofquery[col]<<endl;
        }
    }
    pdist = (float*)mxGetPr(dist);
    double complexity_sum = 0;
    for(size_t r=0; r<ivf.ndocs; r++){
        complexity_sum += ivf.psumofdoc[r];
    }
    for(size_t col=0; col<nqueries; col++){
        double similarity_sum = 0;
        for(size_t r=0; r<ivf.ndocs; r++)
            similarity_sum += pdist[r];
        if(verbose)
            cout<<"complexity "<<complexity_sum<<" similarity "<<similarity_sum<<" ratio "<<complexity_sum/(similarity_sum+0.0001)<<endl;
        float self_dist, min_dist, max_dist;
        min_dist = FLT_MAX; max_dist = FLT_MIN;
        self_dist = psumofquery[col] - psumofquery[col]*w*complexity_sum/(similarity_sum+0.0001);
        vector<size_t> zero_idx;
        for(size_t r=0; r<ivf.ndocs; r++){
            if (pdist[r]>0){ //it has some collsions
                pdist[r] = ivf.psumofdoc[r] - pdist[r]*w*complexity_sum/(similarity_sum+0.0001);
                if (pdist[r]>max_dist)
                    max_dist = pdist[r];
                if (pdist[r]<min_dist)
                    min_dist = pdist[r];
            }
            else
                zero_idx.push_back(r);
        }
        for(size_t t=0; t<zero_idx.size();t++)
            pdist[zero_idx[t]] = max_dist+1;

        //cout<<"orginial"<<endl;
        //cout<<"min_dist:"<<min_dist<<" self_dist:"<<self_dist<<" max_dist:"<<max_dist<<endl;
        assert(min_dist>=self_dist);
        // normalize to [0,2]
        // first make sure the smallest dist > 0
        if(self_dist<0){
            for(size_t r=0; r<ivf.ndocs; r++)
                pdist[r] -= self_dist;
            max_dist -= self_dist;
            self_dist = 0;
        }
        max_dist++;
        //cout<<"computing"<<endl;
        //cout<<"self_dist:"<<self_dist<<" max_dist:"<<max_dist<<endl;
        for(size_t r=0; r<ivf.ndocs; r++)
            pdist[r] = 2*(pdist[r]-self_dist)/max_dist;

        pdist += ivf.ndocs;
    }
    delete []psumofquery;
    return dist;
}

mxArray* ivSearchInvFile_l2autoasym(ivFile& ivf, const mxArray* queries, double w, bool verbose)
{
    if(verbose) cout<<"l2autoasym"<<endl;
    size_t nwords, nqueries;
    nwords = mxGetM(queries);
    nqueries = mxGetN(queries);
    assert(nwords == ivf.nwords);
    mxClassID classID = mxGetClassID(queries);
    assert(classID == mxDOUBLE_CLASS);
    assert(mxIsSparse(queries));
    unsigned int progress = 100;
    unsigned int percentage = nqueries/progress;
    unsigned int verb_thre = 0;
    double avgsqrsumofdoc = 0;
    if(verbose){
        double nz = mxGetNzmax(queries);
        cout<<"query sparsity: "<<nz/nwords/nqueries<<endl;

        for(int i=0; i<ivf.ndocs; i++)
            avgsqrsumofdoc += ivf.psqrsumofdoc[i];

        avgsqrsumofdoc /= ivf.ndocs;
    }
    double* query_pr = (double*) mxGetData(queries);
    size_t* query_jc = (size_t*) mxGetJc(queries);
    size_t* query_ir = (size_t*) mxGetIr(queries);
    double* dist1 = new double[nqueries*ivf.ndocs];
    double* dist2 = new double[nqueries*ivf.ndocs];
    double* pdist1, *pdist2;
    pdist1 = dist1; pdist2 = dist2;
    for (size_t i=0; i<nqueries; i++){
        memcpy(pdist2,ivf.psqrsumofdoc,sizeof(double)*ivf.ndocs);

        //get the staring index for this column
        size_t rstart = query_jc[i];
        size_t rend = query_jc[i+1];
        double sqrsumofquery = 0;
        if (isnormal(*(query_pr+rstart))){
            //column has something
            for (size_t r=rstart; r<rend; r++)
            {
                float tf = float(*(query_pr+r));
                sqrsumofquery += tf*tf;
            }
        }
        for(size_t j=0; j<ivf.ndocs; j++)
            pdist1[j]=sqrsumofquery;

        if(verbose)
            cout<<"avgsqrsumofdoc/sqrsumofquery: "<<avgsqrsumofdoc/sqrsumofquery<<endl;

        pdist1 += ivf.ndocs;
        pdist2 += ivf.ndocs;
    }
    pdist1 = dist1; pdist2 = dist2;
    for(size_t col=0; col<nqueries; col++)
    {
        //get the staring index for this column
        size_t rstart = query_jc[col];
        size_t rend = query_jc[col+1];
        if (isnormal(*(query_pr+rstart))){
            //column has something
            for (size_t r=rstart; r<rend; r++)
            {
                float tf = float(*(query_pr+r));
                size_t word_id = *(query_ir+r);
                for (ivWordDocIt it=ivf.words[word_id].begin(); it!=ivf.words[word_id].end(); it++){
                    float inter = min(tf,it->termFreq);
                    float sqrinter = inter*inter;
                    pdist1[it->docID] += (-2*tf*inter + sqrinter);
                    pdist2[it->docID] += (-2*it->termFreq*inter + sqrinter);
                }
            }
        }
        pdist1 += ivf.ndocs;
        pdist2 += ivf.ndocs;
    }
    mxArray* dist = mxCreateNumericMatrix(ivf.ndocs, nqueries, mxSINGLE_CLASS, mxREAL);
    float *pdist = (float*)mxGetPr(dist);
    pdist1 = dist1; pdist2 = dist2;
    for(size_t col=0; col<nqueries; col++){
        double similarity_sum = 0;
        double complexity_sum = 0;
        for(size_t r=0; r<ivf.ndocs; r++){
            pdist1[r] = sqrt(pdist1[r]);
            pdist2[r] = sqrt(pdist2[r]);
            similarity_sum += pdist1[r];
            complexity_sum += pdist2[r];
        }
        if(verbose)
            cout<<"complexity "<<complexity_sum<<" similarity "<<similarity_sum<<" ratio "<<complexity_sum/(similarity_sum+0.0001)<<endl;
        for(size_t r=0; r<ivf.ndocs; r++)
            pdist[r] = pdist1[r]*w*complexity_sum/(similarity_sum+0.0001)+ pdist2[r];
        pdist += ivf.ndocs;
        pdist1 += ivf.ndocs;
        pdist2 += ivf.ndocs;
    }
    delete []dist1;
    delete []dist2;
    return dist;
}

mxArray* ivSearchInvFile_l2autoasymnew(ivFile& ivf, const mxArray* queries, double w, bool verbose)
{
    if(verbose) cout<<"l2autoasymnew"<<endl;
    size_t nwords, nqueries;
    nwords = mxGetM(queries);
    nqueries = mxGetN(queries);
    assert(nwords == ivf.nwords);
    mxArray* dist = mxCreateNumericMatrix(ivf.ndocs, nqueries, mxSINGLE_CLASS, mxREAL);
    float *pdist = (float*)mxGetPr(dist);
    mxClassID classID = mxGetClassID(queries);
    assert(classID == mxDOUBLE_CLASS);
    assert(mxIsSparse(queries));
    unsigned int progress = 100;
    unsigned int percentage = nqueries/progress;
    unsigned int verb_thre = 0;
    if(verbose){
        double nz = mxGetNzmax(queries);
        cout<<"query sparsity: "<<nz/nwords/nqueries<<endl;
    }
    double* query_pr = (double*) mxGetData(queries);
    size_t* query_jc = (size_t*) mxGetJc(queries);
    size_t* query_ir = (size_t*) mxGetIr(queries);
    for(size_t col=0; col<nqueries; col++)
    {
        //get the staring index for this column
        size_t rstart = query_jc[col];
        size_t rend = query_jc[col+1];
        if (isnormal(*(query_pr+rstart))){
            //column has something
            for (size_t r=rstart; r<rend; r++)
            {
                float tf = float(*(query_pr+r));
                size_t word_id = *(query_ir+r);
                for (ivWordDocIt it=ivf.words[word_id].begin(); it!=ivf.words[word_id].end(); it++){
                    float inter = min(tf,it->termFreq);
                    pdist[it->docID] += inter*inter;
                }
            }
        }
        pdist += ivf.ndocs;
    }
    pdist = (float*)mxGetPr(dist);
    double complexity_sum = 0;
    for(size_t r=0; r<ivf.ndocs; r++){
        complexity_sum += sqrt(ivf.psqrsumofdoc[r]);
    }
    for(size_t col=0; col<nqueries; col++){
        double similarity_sum = 0;
        for(size_t r=0; r<ivf.ndocs; r++)
            similarity_sum += sqrt(pdist[r]);
        if(verbose)
            cout<<"complexity "<<complexity_sum<<" similarity "<<similarity_sum<<" ratio "<<complexity_sum/(similarity_sum+0.0001)<<endl;
        for(size_t r=0; r<ivf.ndocs; r++)
            pdist[r] = sqrt(ivf.psumofdoc[r]) - sqrt(pdist[r])*w*complexity_sum/(similarity_sum+0.0001);
        pdist += ivf.ndocs;
    }
    return dist;
}

mxArray* ivSearchInvFile_l2autoasymnewsqr(ivFile& ivf, const mxArray* queries, double w, bool verbose)
{
    if(verbose) cout<<"l2autoasymnewsqr"<<endl;
    size_t nwords, nqueries;
    nwords = mxGetM(queries);
    nqueries = mxGetN(queries);
    assert(nwords == ivf.nwords);
    mxArray* dist = mxCreateNumericMatrix(ivf.ndocs, nqueries, mxSINGLE_CLASS, mxREAL);
    float *pdist = (float*)mxGetPr(dist);
    mxClassID classID = mxGetClassID(queries);
    assert(classID == mxDOUBLE_CLASS);
    assert(mxIsSparse(queries));
    unsigned int progress = 100;
    unsigned int percentage = nqueries/progress;
    unsigned int verb_thre = 0;
    if(verbose){
        double nz = mxGetNzmax(queries);
        cout<<"query sparsity: "<<nz/nwords/nqueries<<endl;
    }
    double* query_pr = (double*) mxGetData(queries);
    size_t* query_jc = (size_t*) mxGetJc(queries);
    size_t* query_ir = (size_t*) mxGetIr(queries);
    for(size_t col=0; col<nqueries; col++)
    {
        //get the staring index for this column
        size_t rstart = query_jc[col];
        size_t rend = query_jc[col+1];
        if (isnormal(*(query_pr+rstart))){
            //column has something
            for (size_t r=rstart; r<rend; r++)
            {
                float tf = float(*(query_pr+r));
                size_t word_id = *(query_ir+r);
                for (ivWordDocIt it=ivf.words[word_id].begin(); it!=ivf.words[word_id].end(); it++){
                    float inter = min(tf,it->termFreq);
                    pdist[it->docID] += inter*inter;
                }
            }
        }
        pdist += ivf.ndocs;
    }
    pdist = (float*)mxGetPr(dist);
    double complexity_sum = 0;
    for(size_t r=0; r<ivf.ndocs; r++){
        complexity_sum += ivf.psqrsumofdoc[r];
    }
    for(size_t col=0; col<nqueries; col++){
        double similarity_sum = 0;
        for(size_t r=0; r<ivf.ndocs; r++)
            similarity_sum += pdist[r];
        if(verbose)
            cout<<"complexity "<<complexity_sum<<" similarity "<<similarity_sum<<" ratio "<<complexity_sum/(similarity_sum+0.0001)<<endl;
        for(size_t r=0; r<ivf.ndocs; r++)
            pdist[r] = ivf.psumofdoc[r] - pdist[r]*w*complexity_sum/(similarity_sum+0.0001);
        pdist += ivf.ndocs;
    }
    return dist;
}


// HAMMING EMBEDDING functions

void ivBuildInvFile_HE(ivFile* ivf, const mxArray* db, size_t docOffset, bool verbose, const mxArray* words_ini, const mxArray* clip_b)
{
    //if(docOffset>=43)
    //    cout << "doc:" << (docOffset+1)<< "num clips=" << ivf->ndocs << endl;
    ivf->nwords = mxGetM(db);
    ivf->ndocs = mxGetN(db)+docOffset;
    ivf->nHE=mxGetM(clip_b);        //HE - size of the binary signature
    //cout << "doc:" << (docOffset+1)<< "num clips=" << ivf->ndocs << endl;
    //allocate vectors
    ivf->words.resize(ivf->nwords);
    mxClassID classID = mxGetClassID(db);
    assert(classID == mxDOUBLE_CLASS);
    assert(mxIsSparse(db));
    classID = mxGetClassID(words_ini);
    assert(classID == mxINT32_CLASS);
    assert(!mxIsSparse(words_ini));
    classID = mxGetClassID(clip_b);
    assert(classID == mxLOGICAL_CLASS);
    assert(!mxIsSparse(clip_b));
    double* db_pr = (double*) mxGetData(db);
    size_t* db_jc = (size_t*) mxGetJc(db);
    size_t* db_ir = (size_t*) mxGetIr(db);
    bool* b_pr= (bool*) mxGetData(clip_b);   //HE - pointer to clip_b data
    int* w_pr= (int*) mxGetData(words_ini);   //HE - pointer to words_ini data

    unsigned int progress = 100;
    unsigned int percentage = ivf->ndocs/progress;
    unsigned int verb_thre = 0;
    if(verbose){
        double nz = mxGetNzmax(db);
        cout<<"database sparsity: "<<nz/ivf->nwords/ivf->ndocs<<endl;
    }
    ivf->psumofdoc = new double[ivf->ndocs];
    ivf->psqrsumofdoc = new double[ivf->ndocs];
    //if(docOffset>=43)
    //    cout << "doc:" << (docOffset+1)<< "num clips=" << ivf->ndocs << endl;
    for(unsigned int col=0; col<mxGetN(db); col++) {
        //get the staring index for this column
        size_t rstart = db_jc[col];
        size_t rend = db_jc[col+1];
        double& sumofdoc = ivf->psumofdoc[col+docOffset];
        sumofdoc = 0;
        double& sqrsumofdoc = ivf->psqrsumofdoc[col+docOffset];
        sqrsumofdoc = 0;
        if (isnormal(*(db_pr+rstart))){
            //column has something
            for (size_t r=rstart; r<rend; r++) {        //each word present in the video
                //if((docOffset>=43)&&((r==4728)))
                //    cout << "r="<< r<< "/"<< rend;
                double tf = *(db_pr+r);
                //if((docOffset>=43)&&((r==4728)))
                //    cout << "; tf="<< tf;
                int rate = *(w_pr+r-rstart+1)-*(w_pr+r-rstart);
                //if((docOffset>=43)&&((r==4728)))
                //    cout << "; rate="<< rate;
                bool* mat=new bool[rate*ivf->nHE];
                for(int j=0; j<rate; j++){
                    for (int i=0; i<ivf->nHE; i++)
                        mat[i+j*ivf->nHE]=b_pr[i+(*(w_pr+r-rstart)-1+j)*ivf->nHE];     //HE - cropping of clip_b for the speciffic word
                }
                //if((docOffset>=43)&&((r==4728)))
                //    cout << "mat built";
                ivWordDoc wordDoc;
                sumofdoc += tf;
                //if((docOffset>=43)&&((r==4728)))
                //    cout << "; tfacum="<< sumofdoc;
                sqrsumofdoc += tf*tf;
                //if((docOffset>=43)&&((r==4728)))
                //    cout << "; sqrtbla="<< sqrsumofdoc;
                wordDoc.termFreq = tf;
                //if((docOffset>=43)&&((r==4728)))
                //    cout << "; tf="<< wordDoc.termFreq;
                wordDoc.docID=col+docOffset;
                //if((docOffset>=43)&&((r==4728)))
                //    cout << "; doc="<< wordDoc.docID;
                wordDoc.HEsignature=new bool[rate*ivf->nHE];
                //if((docOffset>=43)&&((r==4728)))
                //    cout << "allocate memory"<< r<< "/"<< rend;
                for(int i=0;i<(rate*ivf->nHE);i++)
                    wordDoc.HEsignature[i]=mat[i];                //HE - signature for each point quantized to that word
                wordDoc.Wrate=rate;                     //HE - unwheighted bow value
                //if((docOffset>=43)&&((r==4728)))
                //    cout << "; rate="<< wordDoc.Wrate;
                ivf->words[*(db_ir+r)].push_back(wordDoc);
 /*             if (r==rstart){                             //debbuging print
                        cout<<"videolut "<<wordDoc.docID<< "; word "<<*(db_ir+r)<<"; val "<<rate<<": HE "<<endl;
                        for(int j=0; j<rate; j++){
                            for (int i=0; i<ivf->nHE; i++)
                                cout<< mat[i+j*ivf->nHE];
                            cout<<endl;
                        }
                }
                if((docOffset>=43)&&((r>rstart)))
                    cout << " doc id:" << wordDoc.docID <<"; word id:" << *(db_ir+r) <<"r="<< r <<endl;
                if((docOffset>=43)&&((r==rstart)||(r==rend-1)))
                    cout << " doc id:" << wordDoc.docID <<"; word id:" << *(db_ir+r) <<"r="<< r <<  "; num words:" << wordDoc.Wrate << "; HE:" << ivf->nHE << endl;

                if((docOffset>=43)&&((r==rstart)||(r==rend-1)))
                    cout<<"mat deleted"<<endl;
                cout<<"word "<<*(db_ir+r)<<": doc "<<col<<" val "<<*(db_pr+r)<<endl;*/
                delete []mat;
            }
        }
        if (verbose){
            if(percentage==0){
                mexPrintf("\r%d/%d,sumofdoc:%.2f",col+1,ivf->ndocs,ivf->psumofdoc[col]);
                mexEvalString("drawnow;"); // to print string immediately.
               // cout<<col+1<<"/"<<ivf->ndocs<<endl;
            }
            else if(col>=verb_thre){
                mexPrintf("\r%d%%,sumofdoc:%.2f",col/percentage,ivf->psumofdoc[col]);
                mexEvalString("drawnow;"); // to print string immediately.
                //cout<<col/percentage*progress<<"%"<<endl;
                verb_thre += percentage;
            }
        }
    }
}


mxArray* ivSearchInvFile_HE_l1(ivFile& ivf, const mxArray* queries, bool verbose, const mxArray* words_ini, const mxArray* query_b, const int th)
{
    size_t nwords, nqueries, nHE;
    nwords = mxGetM(queries);
    nqueries = mxGetN(queries);
    nHE=mxGetM(query_b);            //HE - size of the binary signature

    assert(nwords == ivf.nwords);
    assert(nHE == ivf.nHE);         //HE - the binary signatures have to have the same dimension
    mxArray* sim = mxCreateNumericMatrix(ivf.ndocs, nqueries, mxSINGLE_CLASS, mxREAL);
    float *psim = (float*)mxGetPr(sim);
    //for(int i=0;i<2098;i=i+200)
    //    cout<<"pointer: "<<(psim+i)<<"; value: "<<psim[i]<<endl;
    //cout<<"pointer: "<<(psim+2097)<<"; value: "<<psim[2097]<<"; queries: "<<nqueries<<endl;
    mxClassID classID = mxGetClassID(queries);
    assert(classID == mxDOUBLE_CLASS);
    assert(mxIsSparse(queries));
    classID = mxGetClassID(words_ini);
    assert(classID == mxINT32_CLASS);
    assert(!mxIsSparse(words_ini));
    classID = mxGetClassID(query_b);
    assert(classID == mxLOGICAL_CLASS);
    assert(!mxIsSparse(query_b));
    unsigned int progress = 100;
    unsigned int percentage = nqueries/progress;
    unsigned int verb_thre = 0;
    if(verbose){
        double nz = mxGetNzmax(queries);
        cout<<"query sparsity: "<<nz/nwords/nqueries<<endl;
    }
    double* query_pr = (double*) mxGetData(queries);
    size_t* query_jc = (size_t*) mxGetJc(queries);
    size_t* query_ir = (size_t*) mxGetIr(queries);
    bool* b_pr= (bool*) mxGetData(query_b);          //HE - pointer to query_b data
    int* w_pr= (int*) mxGetData(words_ini);   //HE - pointer to words_ini data
    for(size_t col=0; col<nqueries; col++)              //each query. In this case (HE), they're provided one by one, so this runs only one time.
    {
        //get the staring index for this column
        size_t rstart = query_jc[col];
        size_t rend = query_jc[col+1];
        if (isnormal(*(query_pr+rstart))){
            //column has something
            for (size_t r=rstart; r<rend; r++)          //each word
            {
                float tf = float(*(query_pr+r));
                int rate = *(w_pr+r-rstart+1)-*(w_pr+r-rstart);
                bool* mat=new bool[rate*nHE];
                for(int j=0; j<rate; j++){
                    for (int i=0; i<nHE; i++)
                        mat[i+j*nHE]=b_pr[i+(*(w_pr+r-rstart)-1+j)*nHE];     //HE - cropping of query_b for the speciffic word
                }
                size_t word_id = *(query_ir+r);
                int mr,om;
                float psimil, oldsim;
                for (ivWordDocIt it=ivf.words[word_id].begin(); it!=ivf.words[word_id].end(); it++){  //distance of that word for each video
                    /*debugging print
                     if (((r==rstart)||(r==rstart+1))&&((it==(ivf.words[word_id].begin()))||(it==(ivf.words[word_id].end()-1))))
                                cout<<"word "<<word_id<<"; videolut "<<it->docID <<"; tf query "<<tf<<" tf video "<<(it->termFreq)<<endl;*/
                    int* dist=new int[rate*it->Wrate];
                    int count=0;
                    int countq=0;
                    int countv=0;
                    for(int i=0;i<rate;i++){              //HE - HE distance between all points query-video for that word.
                        for(int j=0;j<it->Wrate;j++){    //look for matches between query vs all points in video
                            int d=0;
                            for(int n=0;n<nHE;n++)
                                d+=mat[n+i*nHE]^it->HEsignature[n+j*nHE];
                            dist[j+i*it->Wrate]=d;
                            //debugging print
                           /*if ((r==rend-1)&&(it==ivf.words[word_id].begin())&&(j==it->Wrate-1)&&(i==rate-1))
                                cout<<"word "<<word_id<<"; videolut "<<it->docID <<";  "<<rate<<" q vs "<<it->Wrate<<"; dist: "<<dist[j+i*it->Wrate]<<endl;*/

                            if (dist[j+i*it->Wrate]<th){             //if HE<predefined th, then there is a match.
                                dist[j+i*it->Wrate]=1;
                                count=1;
                            }else
                                dist[j+i*it->Wrate]=0;
                        }
                        countq+=count;                  //num of matched points in the query.
                        count=0;
                    }
                    for(int i=0;i<it->Wrate;i++){    //look for matches between each point in the video
                        for(int j=0;j<rate;j++){              //and all points in query
                            if (dist[j+i*it->Wrate]==1)          //there is at least one match.
                                count=1;
                           }
                        countv+=count;                  //num of matched points in the query.
                        count=0;
                    }


                    psimil = min(countq,countv)*min(it->termFreq,tf)/min(rate,it->Wrate);
                    psim[it->docID]+=psimil;
                    /*oldsim=min(it->termFreq,tf);
                    mr=min(countq,countv);
                    om=min(rate,it->Wrate);
                    if (((r==rstart)||(r==rstart+1))&&((it==(ivf.words[word_id].begin()))||(it==(ivf.words[word_id].end()-1)))){                             //debbuging print
                        cout<<"word "<<word_id<<"; videolut "<<it->docID <<"; val sim: "<<psimil <<"; old sim"<< oldsim<<"; HE match"<< mr<<"; old match"<< om<<"; sim total: "<<psim[it->docID]<<": "<<rate<<" q vs "<<it->Wrate<<" -> "<<countq<<" to "<<countv<<endl;
                        for(int j=0; j<rate; j++){
                            for (int i=0; i<it->Wrate; i++)
                                cout<< dist[i+j*it->Wrate];
                            cout<<endl;
                        }
                    }*/

                    delete []dist;
                }
            }
        }
 /*       for(int i=0;i<2098;i=i+16){
            if(i%6==0) cout<<endl;
            cout<<"; ID:"<<i<<"; value: "<<psim[i];
        }
        cout<<"; value 2097:"<<psim[2097]<<endl;*/
        if (nqueries>1) psim += ivf.ndocs;
        if (verbose && nqueries>1){
            if(percentage==0){
                mexPrintf("\r%d/%d",col+1,nqueries);
                mexEvalString("drawnow;"); // to print string immediately.
                //cout<<col+1<<"/"<<nqueries<<endl;
            }
            else if(col>=verb_thre){
                mexPrintf("\r%d%%",col/percentage);
                mexEvalString("drawnow;"); // to print string immediately.
                //cout<<col/percentage<<"%"<<endl;
                verb_thre += percentage;
            }
        }
    }
   // cout<<"ok query"<<nqueries<<" for num videos: "<<ivf.ndocs<<endl;
    return sim;
}
