/*   asl_functions.h various functions for the manipulation of ASL data

      Michael Chappell - FMIRB Image Analysis Group

      Copyright (C) 2009 University of Oxford */

/*   CCOPYRIGHT   */

#if !defined(asl_functions_h)
#define asl_functions_h

#include <string>
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;

namespace OXASL {
 
  void data2stdform(Matrix& datamtx, vector<Matrix>& asldata, int ntis, bool isblocked, bool ispairs, bool blockpairs);
  void stdform2data(vector<Matrix>& asldata, Matrix& datareturn, bool outblocked, bool outpairs);

  // separate the pairs in the data (into seprate standard form items)
  void separatepairs(vector<Matrix>& asldata, vector<Matrix>& asldataodd, vector<Matrix>& asldataeven, bool blockpairs);
  void mergepairs(vector<Matrix>& asldata, vector<Matrix>& asldataodd, vector<Matrix>&asldataeven);

  // mean the data at each TI
  void timeans(vector<Matrix>& asldata, vector<Matrix>& aslmean);
  // output data having taken mean at each TI
  void timeanout(vector<Matrix>& asldata,  volume<float>& mask, string fname, bool outpairs);

  // output data split into individual TIs
  void splitout(vector<Matrix>& asldata, volume<float>& mask, string froot);

  // output epochs of data
  void epochout(vector<Matrix>& asldata, volume<float>& mask, string froot, int epadv, int epol, bool outpairs, bool tiunit=false);
  // generate epochs
  void genepochs(vector<Matrix>& asldata, vector<Matrix>& epochreturn, int epadv, int epol);
  // generate TI epochs
  void gentiepochs(Matrix& asldata, vector<Matrix>& epochreturn, int epadv, int epol);

  //output the result of deconvolution with an AIF
  void deconvout(vector<Matrix>& asldata, volume<float>& mask, Matrix& aif, string fname);
  // do SVD convoloution
  ReturnMatrix SVDdeconv(const Matrix& data, const Matrix& aif);
  // create a (simple) convolution matrix
  ReturnMatrix convmtx(const ColumnVector& invec);
}

#endif
