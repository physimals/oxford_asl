/*   asl_functions.h various functions for the manipulation of ASL data

    Michael Chappell - FMRIB Image Analysis Group

    Moss Zhao - IBME Quantitative Biomedical Inference (QuBIc) Group

    Copyright (C) 2015 University of Oxford  */

/*   CCOPYRIGHT   */

#if !defined(asl_functions_h)
#define asl_functions_h

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "armawrap/newmat.h"
#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"

namespace OXASL
{
void data2stdform(NEWMAT::Matrix &datamtx, std::vector<NEWMAT::Matrix> &asldata, int ntis, std::vector<int> nrpts, bool isblocked, bool ispairs, bool blockpairs);
void stdform2data(std::vector<NEWMAT::Matrix> &asldata, NEWMAT::Matrix &datareturn, bool outblocked, bool outpairs);

// separate the pairs in the data (into seprate standard form items)
void separatepairs(std::vector<NEWMAT::Matrix> &asldata, std::vector<NEWMAT::Matrix> &asldataodd, std::vector<NEWMAT::Matrix> &asldataeven);
void mergepairs(std::vector<NEWMAT::Matrix> &asldata, std::vector<NEWMAT::Matrix> &asldataodd, std::vector<NEWMAT::Matrix> &asldataeven);

// mean the data at each TI
void timeans(std::vector<NEWMAT::Matrix> &asldata, std::vector<NEWMAT::Matrix> &aslmean);
// output data having taken mean at each TI
void timeanout(std::vector<NEWMAT::Matrix> &asldata, NEWIMAGE::volume<float> &mask, std::string fname, bool outpairs);

// output data split into individual TIs
void splitout(std::vector<NEWMAT::Matrix> &asldata, NEWIMAGE::volume<float> &mask, std::string froot);

// output epochs of data
void epochout(std::vector<NEWMAT::Matrix> &asldata, NEWIMAGE::volume<float> &mask, std::string froot, int epadv, int epol, bool outpairs, bool tiunit = false);
// generate epochs
void genepochs(std::vector<NEWMAT::Matrix> &asldata, std::vector<NEWMAT::Matrix> &epochreturn, int epadv, int epol);
// generate TI epochs
void gentiepochs(NEWMAT::Matrix &asldata, std::vector<NEWMAT::Matrix> &epochreturn, int epadv, int epol);

//output the result of deconvolution with an AIF
void deconvout(std::vector<NEWMAT::Matrix> &asldata, NEWIMAGE::volume<float> &mask, NEWMAT::Matrix &aif, std::string fname);
// do SVD convoloution
NEWMAT::ReturnMatrix SVDdeconv(const NEWMAT::Matrix &data, const NEWMAT::Matrix &aif);
// create a (simple) convolution matrix
NEWMAT::ReturnMatrix convmtx(const NEWMAT::ColumnVector &invec);

// Partial volume correction functions
// function to perform partial volume correction by linear regression
void pvcorr_LR(std::vector<NEWMAT::Matrix> &data_in, int ndata_in, NEWIMAGE::volume<float> &mask, NEWIMAGE::volume<float> &pv_map_gm, NEWIMAGE::volume<float> &pv_map_wm, int kernel, std::vector<NEWMAT::Matrix> &data_out, bool outblocked, bool outpairs, std::vector<int> nrpts, bool isblocked, bool ispairs, bool blockpairs);

// PV correction using linear regression (Asllani's method)
NEWIMAGE::volume<float> correct_pv_lr(const NEWIMAGE::volume<float> &data_in, const NEWIMAGE::volume<float> &mask, const NEWIMAGE::volume<float> &pv_map_gm, const NEWIMAGE::volume<float> &pv_map_wm, int kernel);

// Function to correct NaN values
NEWIMAGE::volume<float> correct_NaN(const NEWIMAGE::volume<float> &data_in);

// Extrapolation functions
// function to extrapolate voxels
void extrapolate(std::vector<NEWMAT::Matrix> &data_in, int ndata_in, int ntis, NEWIMAGE::volume<float> &mask, int neighbour_size, std::vector<NEWMAT::Matrix> &data_out, bool outblocked, bool outpairs, std::vector<int> nrpts, bool isblocked, bool ispairs, bool blockpairs);

// Function to do spiral search on 2D images and extrapolate edge voxels
NEWMAT::Matrix extrapolate_avg(NEWMAT::Matrix data_in, NEWMAT::Matrix mask, int neighbour_size);
}

#endif
