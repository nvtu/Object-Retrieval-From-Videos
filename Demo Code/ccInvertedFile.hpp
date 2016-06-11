#ifndef CC_INVERTEDFILE
#define CC_INVERTEDFILE

#include "mex.h"

#include <iostream>
#include <string>
#include <utility>
#include <deque>
#include <stdint.h>
using namespace std;

typedef struct _ivWordDoc{
    unsigned int  docID;
    float termFreq;
//HE parameters:
    int Wrate;
    bool* HEsignature;
} ivWordDoc;
//typedef pair<uint,float> ivWordDoc;
typedef deque<ivWordDoc> ivWord;
typedef deque<ivWordDoc>::iterator ivWordDocIt;
typedef deque<ivWord> ivWordList;
typedef ivWordList::iterator ivWordIt;

//------------------------------------------------------------------------
//Inverted file structure

class ivFile
{
public:
  //array of word entries
  ivWordList words;
  unsigned int nwords;
  unsigned int ndocs;
  unsigned int nHE;
  ivFile()
  {
      nwords = 0;
      ndocs = 0;
      psumofdoc = NULL;
      psqrsumofdoc = NULL;
  }
  ~ivFile()
  {
      if (psumofdoc!=0){
          delete []psumofdoc;
          psumofdoc = NULL;
      }
      if (psqrsumofdoc!=0){
          delete []psqrsumofdoc;
          psqrsumofdoc = NULL;
      }
      int t,i;
      int s = words.size()/10;
      t=i=0;
      for( ivWordIt wit=words.begin(); wit!=words.end(); wit++)
      {
//          if (i++>=t){
//              mexPrintf("\r%d\%",i/s*10);
//              mexEvalString("drawnow;"); // to print string immediately.
//              t+=s;
//          }
          ivWord().swap(*wit);
      }
      ivWordList().swap(words);
  }
  double *psumofdoc;
  double *psqrsumofdoc;
  double sparsity;
//  //write to file
//  void save(string filename);
//
//  //load from file
//  void load(string filename);
//
//  void display();
//
//  //clears the memory
//  void clear();

  //fill the inverted file with input counts
  //
  // data     - the input data, with one data vector per input consisting of
  //            all the word labels for its tokens
  friend void ivBuildInvFile(ivFile* ivf, const mxArray* db, size_t docOffset=0, bool verbose=true);
  friend mxArray* ivComputeSelfSim(ivFile& ivf, const mxArray* db,  bool verbose=true);
  friend mxArray* ivComputeSelfSim_l1(ivFile& ivf, const mxArray* db,  bool verbose=true);
  friend mxArray* ivComputeSelfSim_l2(ivFile& ivf, const mxArray* db,  bool verbose=true);
  friend void ivBuildInvFile_HE(ivFile* ivf, const mxArray* db, size_t docOffset=0, bool verbose=true, const mxArray* words_ini=NULL, const mxArray* clip_b=NULL);


  // search the inverted file for the closest document
  //
  // data     - the input data, with one vector per iput
  // overlapOnly - return only those documents with overlapping words
  // k        - no. of output required per document, if 0 then return everything
  // scorelists   - a list of ivNOdeList to hold the results
  friend mxArray* ivSearchInvFile_l1(ivFile& ivf, const mxArray* queries, bool verbose=true);
  friend mxArray* ivSearchInvFile_l2(ivFile& ivf, const mxArray* queries, bool verbose=true);
  friend mxArray* ivSearchInvFile_bridge_l1(ivFile& ivf, const mxArray* queries, const mxArray* query_bridges, bool verbose=true);
  friend mxArray* ivSearchInvFile_bridge_l2(ivFile& ivf, const mxArray* queries, const mxArray* query_bridges, bool verbose=true);
  friend mxArray* ivSearchInvFile_asym(ivFile& ivf, const mxArray* queries, double w1, double w2, bool verbose=true);
  friend mxArray* ivSearchInvFile_autoasym(ivFile& ivf, const mxArray* queries, double w, bool verbose=true);
  friend mxArray* ivSearchInvFile_l2asym(ivFile& ivf, const mxArray* queries, double w1, double w2, bool verbose=true);
  friend mxArray* ivSearchInvFile_l2autoasym(ivFile& ivf, const mxArray* queries, double w, bool verbose=true);
  friend mxArray* ivSearchInvFile_l2autoasymnew(ivFile& ivf, const mxArray* queries, double w, bool verbose=true);
  friend mxArray* ivSearchInvFile_l2autoasymnewsqr(ivFile& ivf, const mxArray* queries, double w, bool verbose=true);
  friend mxArray* ivSearchInvFile_HE_l1(ivFile& ivf, const mxArray* queries, bool verbose=true, const mxArray* words_ini=NULL, const mxArray* query_b=NULL, const int th=22);

  friend void ivCleanInvFile(ivFile* ivf){delete ivf;}
};


#endif
