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
#define dbIn      prhs[1]
#define verbIn    prhs[2]
#define normIn    prhs[3]
//outputs          
#define simOut    plhs[0]
          
         

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  if (nrhs!=4 && nrhs!=3)	mexErrMsgTxt("3-4 inputs required");

  //make the ivFile object  
  ivFile* ivfile = *(ivFile**)mxGetData(ivFileIn);
  //check if not passed in an object, then create a new one
  assert (ivfile != NULL);

  bool verbose = *(bool*)mxGetData(verbIn);
  mxArray* sim;
  if (nrhs==3)
      sim = ivComputeSelfSim(*ivfile, dbIn, verbose);
  else{
      double *pPara = (double*)mxGetPr(normIn);
      int opt = round(pPara[0]);
      switch (opt){
          case 1:
              sim = ivComputeSelfSim_l1(*ivfile, dbIn, verbose);
              break;
          case 2:
              sim = ivComputeSelfSim_l2(*ivfile, dbIn, verbose);
              break;
      } 
  }
  simOut = sim;
}

