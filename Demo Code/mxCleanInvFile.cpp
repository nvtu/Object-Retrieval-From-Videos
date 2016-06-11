// Author: Mohamed Aly <malaa at vision d0t caltech d0t edu>
// Date: October 6, 2010

#include "mex.h"

#include <string>
#include <cassert>
#include <cstring>
#include "ccInvertedFile.hpp"

//inputs
#define ivFileIn  prhs[0]
          
         

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
   
  if (nrhs !=1)	mexErrMsgTxt("1 inputs required");

  //make the ivFile object  
  ivFile* ivfile = *(ivFile**)mxGetData(ivFileIn);
  //check if not passed in an oVbject, then create a new one
  if (ivfile != NULL)
    ivCleanInvFile(ivfile);
}

