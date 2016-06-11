// Author: Mohamed Aly <malaa at vision d0t caltech d0t edu>
// Date: October 6, 2010

#include "mex.h"

#include <string>
#include <cassert>
#include <cstring>
#include <climits>
#include "ccInvertedFile.hpp"

//inputs
#define ivFileIn  prhs[0]
#define dataIn    prhs[1]
#define docOffIn  prhs[2]
#define verbIn    prhs[3]
#define wordsIn   prhs[4]
#define HEIn      prhs[5]

//outputs          
#define ivFileOut plhs[0]
                   

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
   
  if (nrhs<2 || nrhs>6)	mexErrMsgTxt("2~6 inputs required");

  //make the ivFile object  
  ivFile* ivfile = *(ivFile**)mxGetData(ivFileIn);
  //check if not passed in an object, then create a new one
  if (ivfile == NULL)
    ivfile = new ivFile;
   

  if (nrhs==2)
      ivBuildInvFile(ivfile, dataIn);
  else{
      size_t docOffset = (size_t) *mxGetPr(docOffIn);
      bool verbose = true;
      if (nrhs>3) verbose = *(bool*)mxGetData(verbIn);
      if(nrhs==6){
      //  cout<<"ivBuildInvFile_HE(ivfile, dataIn, docOffset, verbose,wordsIn,HEIn);"<<endl;
        ivBuildInvFile_HE(ivfile, dataIn, docOffset, verbose,wordsIn,HEIn);
      }else
       //           cout<<"ivBuildInvFile_HE(ivfile, dataIn, docOffset, verbose);"<<endl;
        ivBuildInvFile(ivfile, dataIn, docOffset, verbose);
  }
      

  //return pointer to ivFile
  mxArray* ret = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
  *(ivFile**) mxGetPr(ret) = ivfile;
  ivFileOut = ret;
}

