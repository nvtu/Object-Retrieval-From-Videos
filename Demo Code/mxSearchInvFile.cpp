// Author: Mohamed Aly <malaa at vision d0t caltech d0t edu>
// Date: October 6, 2010

#include "mex.h"
#include <cmath>
#include <string>
#include <cassert>
#include <cstring>
#include <climits>
#include "ccInvertedFile.hpp"

//inputs
#define ivFileIn  prhs[0]
#define queryIn   prhs[1]
#define paraIn    prhs[2]
#define verbIn    prhs[3]
#define bridges   prhs[4]
#define wordsIn   prhs[4]
#define HEIn      prhs[5]
#define thIn      prhs[6]


//outputs          
#define simOut    plhs[0]
          
         

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
   
  if ((nrhs!=4)&&(nrhs!=5)&&(nrhs!=7))	mexErrMsgTxt("4, 5 or 7 inputs required");

  //make the ivFile object  
  ivFile* ivfile = *(ivFile**)mxGetData(ivFileIn);
  //check if not passed in an object, then create a new one
  assert (ivfile != NULL);

  double *pPara = (double*)mxGetPr(paraIn);
  int opt = round(pPara[0]);
  bool verbose = *(bool*)mxGetData(verbIn);
  mxArray* sim;
  switch (opt){
      case 1: 
          if(nrhs==4)
            sim = ivSearchInvFile_l1(*ivfile, queryIn, verbose);
          else if(nrhs==7){
          	int th = *(int*)mxGetData(thIn);
            sim = ivSearchInvFile_HE_l1(*ivfile, queryIn, verbose, wordsIn, HEIn, th);
          }
          break;
      case 2: 
          sim = ivSearchInvFile_l2(*ivfile, queryIn, verbose);
          break;
      case 3: 
          sim = ivSearchInvFile_asym(*ivfile, queryIn, pPara[1], pPara[2], verbose);
          break;
      case 4: 
          sim = ivSearchInvFile_autoasym(*ivfile, queryIn, pPara[1], verbose);
          break;
      case 5: 
          sim = ivSearchInvFile_l2asym(*ivfile, queryIn, pPara[1], pPara[2], verbose);
          break;
      case 6: 
          sim = ivSearchInvFile_l2autoasym(*ivfile, queryIn, pPara[1], verbose);
          break;
      case 7: 
          sim = ivSearchInvFile_l2autoasymnew(*ivfile, queryIn, pPara[1], verbose);
          break;
      case 8: 
          sim = ivSearchInvFile_l2autoasymnewsqr(*ivfile, queryIn, pPara[1], verbose);
          break;
      case 9: 
          sim = ivSearchInvFile_bridge_l1(*ivfile, queryIn, bridges, verbose);
          break;
      case 10: 
          sim = ivSearchInvFile_bridge_l2(*ivfile, queryIn, bridges, verbose);
          break;
  }

  simOut = sim;
}

