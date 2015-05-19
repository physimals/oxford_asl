/*   asl_functions.cc various functions for the manipulation of ASL data

    Michael Chappell - FMRIB Image Analysis Group

    Moss Zhao - IBME Quantitative Biomedical Inference (QuBIc) Group

    Copyright (C) 2015 University of Oxford  */

/*   CCOPYRIGHT   */

#include "asl_functions.h"

namespace OXASL {

  void data2stdform(Matrix& datamtx, vector<Matrix>& asldata, int ntis, bool isblocked, bool ispairs) {
    int nvox=datamtx.Ncols();
    int nmeas=datamtx.Nrows()/ntis;
    int nrpts;
    if (ispairs) nrpts=nmeas/2;
    else nrpts=nmeas;

    //cout << nmeas << " " << nrpts << " " << endl;

    if (isblocked)
      {
	Matrix thisti(nmeas,nvox);
	for (int ti=1; ti<=ntis; ti++) 
	  {
	    thisti=0;
	    //asldata[ti-1].ReSize(nvox,nmeas);
	    //extract the measurements for this TI
	    for (int i=1; i<=nrpts; i++)
	      {
		if (ispairs) {
		  thisti.Row(2*i-1) = datamtx.Row(2*(i-1)*ntis+2*ti-1);
		  thisti.Row(2*i) = datamtx.Row(2*(i-1)*ntis+2*ti);
		}
		else {
		  thisti.Row(i) = datamtx.Row((i-1)*ntis+ti);
		}
	      }
	  
	    asldata.push_back(thisti);
	  }

	/*
	//deal with orphan data at end (if present)
	if (datamtx.Nrows() > nmeas*ntis) {
	  int i=nrpts*2*ntis+1;
	  int ti=0;
	  while (i<= datamtx.Nrows()) {
	    if (ispairs) {
	      asldata[ti] = asldata[ti] & datamtx.Row(i);
	      i++;
	      if (i > datamtx.Nrows()) throw Exception("Cannot process this dataset as pairs: odd number of spare measurements at end");
	      asldata[ti] = asldata[ti] & datamtx.Row(i);
	    }
	    else {
	      asldata[ti] = asldata[ti] & datamtx.Row(i);
	    }
	    i++;
	    ti++;
	  }
	}
	*/
	
      }
    else 
      {

	for (int ti=1; ti<=ntis; ti++) 
	  //asldata[ti-1].ReSize(nvox,nmeas);
	  asldata.push_back(datamtx.Rows((ti-1)*nmeas+1,ti*nmeas));

	if (datamtx.Nrows() > nmeas*ntis) throw Exception("Orphaned data found at end of file - this is not logical when data is in TI blocks");
      }
  }

  void stdform2data(vector<Matrix>& asldata, Matrix& datareturn, bool outblocked, bool outpairs) {
    int ntis = asldata.size();
    int nvox = asldata[0].Ncols();
    int nmeas = asldata.back().Nrows(); //safer to determine this from the very last TI (in case of orphan measurements when nodiscard is turned on)
    int ninc=1;
    if (outpairs) ninc=2;
    
    if (outblocked) {
      datareturn.ReSize(ntis*nmeas,nvox);
      int idx=1;
      for (int m=1; m<=nmeas; m+=ninc)
	for (int ti=1; ti<=ntis; ti++) {
	  datareturn.Row(idx) = asldata[ti-1].Row(m);
	  idx++;
	  if (outpairs) {
	    datareturn.Row(idx) = asldata[ti-1].Row(m+1);
	    idx++;
	  }
	}
      assert(idx-1==ntis*nmeas);

      /*
      //deal with orphans
      if (asldata.front().Nrows() > nmeas) {
	int ti=0;
	int finalm = nmeas+ninc;
	while (asldata[ti].Nrows() > nmeas) {
	  datareturn.Row(idx) = asldata[ti].Row(finalm);
	  idx++;
	  if (outpairs) {
	    datareturn.Row(idx) = asldata[ti].Row(finalm+1);
	    idx++;
	  }
	  ti++;
	}
      }
      */
    }
    else {
      datareturn = asldata[0];
      for (int ti=1; ti<ntis; ti++) {
	datareturn &= asldata[ti];
      }
    }
  }

  void separatepairs(vector<Matrix>& asldata, vector<Matrix>& asldataodd, vector<Matrix>& asldataeven) {
    int ntis = asldata.size();
    int nmeas = asldata[0].Nrows();
    int nrpts=nmeas/2; //if we are using this function then the data must contain pairs

    // just in case the vectors are not empty to start with
    asldataodd.clear();
    asldataeven.clear();

    int idx;

    for (int ti=0; ti<ntis; ti++) {

	Matrix oddmtx;
	oddmtx = asldata[ti].Row(1);
	Matrix evenmtx;
	evenmtx = asldata[ti].Row(2);

	for (int r=2; r<=nrpts; r++) {
	  idx=(r-1)*2+1;
	  oddmtx &= asldata[ti].Row(idx);
	  evenmtx &= asldata[ti].Row(idx+1);
	}
	asldataodd.push_back(oddmtx);
	asldataeven.push_back(evenmtx);
      }
  }

  void mergepairs(vector<Matrix>& asldata, vector<Matrix>& asldataodd, vector<Matrix>&asldataeven) {
    int ntis = asldataodd.size();
    int nmeas = asldataodd[0].Nrows();
    int nrpts=nmeas; //asldataodd does not contain pairs

    asldata.clear(); //make sure this is clear

    for (int ti=0; ti<ntis; ti++) {
      Matrix aslmtx;
      aslmtx = asldataodd[ti].Row(1);
      aslmtx &= asldataeven[ti].Row(1);
      for (int r=2; r<=nrpts; r++) {
	aslmtx &= asldataodd[ti].Row(r);
	aslmtx &= asldataeven[ti].Row(r);
      }
      asldata.push_back(aslmtx);
    }
  }

  void timeans(vector<Matrix>& asldata, vector<Matrix>& meanreturn) {
    int ntis = asldata.size();
    int nvox = asldata[0].Ncols();
    Matrix meanti(ntis,nvox);
      for (int n=0; n<ntis; n++) 
	{
	  meanreturn.push_back(mean(asldata[n],1));
	}
  }

  void timeanout(vector<Matrix>& asldata, volume<float>& mask, string fname, bool outpairs) {
      cout << "Outputting ASL data mean at each TI" << endl;
      
	vector<Matrix> meanti;
      if (!outpairs) {
	timeans(asldata,meanti);
      }
      else {
	//need to preserve pairs in the data - separate
	vector<Matrix> asldataodd;
	vector<Matrix> asldataeven;
	separatepairs(asldata,asldataodd,asldataeven);

	//take mean of odds and evens
	vector<Matrix> meanodd;
	timeans(asldataodd,meanodd);
	vector<Matrix> meaneven;
	timeans(asldataeven,meaneven);

	//recombine
	mergepairs(meanti,meanodd,meaneven);
      }

      Matrix meantimat;
      stdform2data(meanti, meantimat, false, false);

      volume4D<float> meanout;
      meanout.setmatrix(meantimat,mask);
      save_volume4D(meanout,fname);
      
    }

  void splitout(vector<Matrix>& asldata, volume<float>& mask, string froot) {
      cout << "Splitting ASL data into files for each TI" << endl;
      int ntis = asldata.size();
      volume4D<float> blockout;
      for (int n=0; n<ntis; n++) 
	{
	  char cstr [5];
	  if (n<10) sprintf(cstr,"00%d",n);
	  else if (n<100) sprintf(cstr,"0%d",n);
	  else if (n<1000) sprintf(cstr,"%d",n);
	  else throw Exception("More than 1000 measurements in this ASL data file, sorry cannot handle this operation");
	  string tino(cstr);
	  
	  blockout.setmatrix(asldata[n],mask);
	  save_volume4D(blockout,froot+tino);
	}
    }

  void genepochs(vector<Matrix>& asldata, vector<Matrix>& epochreturn, int eplen, int epol) {
      
      int ntis = asldata.size();
      int nvox = asldata[0].Ncols();
      int nmeas = asldata[0].Nrows();

      epochreturn.clear();

      if (epol>=eplen) throw Exception("The epoch overlap may not exceed or equal the length of the epoch");
      int epadv = eplen-epol;

      int e=1;
      Matrix epoch_temp(ntis,nvox);
      while ((e-1)*epadv+eplen<=nmeas)
	{
	  //Matrix ti_temp(eplen,nvox);
	  //go through the TIs
	  for (int ti=1; ti<=ntis; ti++) {
	    epoch_temp.Row(ti) = mean(asldata[ti-1].Rows((e-1)*epadv+1,(e-1)*epadv+eplen),1);
	  }

	  epochreturn.push_back(epoch_temp);

	  e++;
	}

      e--; //correct for final e++
      int unused=nmeas-( (e-1)*epadv+eplen );
      if (unused>0) cout << "Number of measurements from end of data discarded: " << unused << endl;
    }

  void gentiepochs(Matrix& asldata, vector<Matrix>& epochreturn, int eplen, int epol) {
      // expects data in matrix form (blocks of repeats)

      int nvox = asldata.Ncols();
      int nmeas = asldata.Nrows();

      epochreturn.clear();

      if (epol>=eplen) throw Exception("The epoch overlap may not exceed or equal the length of the epoch");
      int epadv = eplen-epol;

      int e=1;
      Matrix epoch_temp(eplen,nvox);
      while ((e-1)*epadv+eplen<=nmeas)
	{
	  epoch_temp = asldata.Rows( (e-1)*epadv+1,(e-1)*epadv+eplen );
	  epochreturn.push_back(epoch_temp);
	  e++;
	}

      e--;
      
      int unused=nmeas-( (e-1)*epadv+eplen );
      if (unused>0) cout << "Number of measurements from end of data discarded: " << unused << endl;
    }


    void epochout(vector<Matrix>& asldata, volume<float>& mask, string froot, int eplen, int epol, bool outpairs, bool tiunit) {
      //Save out epochs of data.
      // tiunit indicates that we want epochs to be done over the Tis rather than over the repeats
    cout << "Ouput ASL data epochs" << endl;
 
    //generate the epochs
    vector<Matrix> theepochs;
    Matrix alldata;
    if (!outpairs) {
      if (!tiunit) {
	genepochs(asldata, theepochs, eplen, epol);
      }
      else {
	// epochs over the TIs need to convert the data to matrix form
	stdform2data(asldata,alldata,true,false);
	gentiepochs(alldata,theepochs,eplen,epol);
      }
    }
    else {
      //need to preserve pairs in the data - separate
      vector<Matrix> asldataodd;
      vector<Matrix> asldataeven;
      separatepairs(asldata,asldataodd,asldataeven);

      vector<Matrix> oddepochs;
      if (!tiunit) {
	genepochs(asldataodd, oddepochs, eplen, epol);
      }
      else {
	stdform2data(asldataodd,alldata,true,false);
	gentiepochs(alldata, oddepochs, eplen, epol);
      }

      vector<Matrix> evenepochs;
      if (!tiunit) {
	genepochs(asldataeven, evenepochs, eplen, epol);
      }
      else {
	stdform2data(asldataodd,alldata,true,false);
	gentiepochs(alldata, evenepochs, eplen, epol);
	
      }

      mergepairs(theepochs,oddepochs,evenepochs);
    }

    int nepoch = theepochs.size();
    volume4D<float> epochout;
    for (int e=0; e<nepoch; e++) {
      char cstr [5];
	  if (e<10) sprintf(cstr,"00%d",e);
	  else if (e<100) sprintf(cstr,"0%d",e);
	  else if (e<1000) sprintf(cstr,"%d",e);
	  else throw Exception("More than 1000 epochs in this ASL data file, sorry cannot handle this operation");
	  string epno(cstr);

	  epochout.setmatrix(theepochs[e],mask);
	  save_volume4D(epochout,froot+epno);
    }
  }

  ReturnMatrix SVDdeconv(const Matrix& data, const Matrix& aif) {
    // do a singular value deconvolution of the data to get residue function

    int nti = data.Nrows();
    int nvox = data.Ncols();
    float truncfac = 0.2;

    // voxelwise SVD deconvolution
    Matrix aifconv; Matrix residue(nti,nvox);
    DiagonalMatrix S; DiagonalMatrix D; Matrix U; Matrix V;
    for (int v=1; v<=nvox; v++) {
      //make convolution matrix
      aifconv = convmtx(aif.Column(v));
      //SVD
      SVD(aifconv,S,U,V);
      // invert the singular values
      D = S.i();
      // truncate (zero all singular values below threshold)
      for (int i=2; i<=D.Nrows(); i++) {
	if (S(i,i) < truncfac*S(1,1)) D(i,i)=0;
      }
      // calculate resdiue
      residue.Column(v) = V*D*U.t()*data.Column(v);
    }
      
    return residue; 
  }

  ReturnMatrix convmtx(const ColumnVector& invec){
    // create a (simple) convolution matrix

    int nentry = invec.Nrows();
    Matrix cmat(nentry,nentry);
    cmat=0.0;
    for (int i=1; i<=nentry; i++) {
      cmat.SubMatrix(i,i,1,i) = ((invec.Rows(1,i)).Reverse()).AsRow();
    }

    return cmat;
  }

  void deconvout(vector<Matrix>& asldata, volume<float>& mask, Matrix& aif, string fname) {
    //perform deconvolution and output the magnitude and residue function

    //take the mean in each TI (we dont want mulitple repeats here)
    vector<Matrix> meandata;
    timeans(asldata, meandata);
    Matrix data;
    stdform2data(meandata,data,false,false); //way to get mean result from standard form into a matrix that can now be processed

    //do the deconvolution
    Matrix resid = SVDdeconv(data,aif);
    // extract magntiude and residue separately
    int nvox = data.Ncols();
    ColumnVector mag(nvox);
    for (int v=1; v<=nvox; v++) {
      mag(v) = (resid.Column(v)).Maximum();
      resid.Column(v) /= mag(v);
    }

    //output 
    volume4D<float> residout;
    residout.setmatrix(resid,mask);
    save_volume4D(residout,fname+"_residuals");

    volume4D<float> magout;
    magout.setmatrix(mag.AsMatrix(1,nvox),mask);
    save_volume4D(magout,fname+"_magntiude");
  }

  // function to perform partial volume correction by linear regression
  void pvcorr_LR(const volume4D<float>& data_in, int ndata_in, const volume<float>& mask, const volume<float>& pv_map, int kernel, volume4D<float>& data_pvcorr) {

    // Clone input data to pv corrected data
    data_pvcorr = data_in;

    // Correct NaN and INF numbers of input mask and pvmap
    volume<float> mask_in_corr(mask.xsize(), mask.ysize(), mask.zsize());
    volume<float> pv_map_in_corr(pv_map.xsize(), pv_map.ysize(), pv_map.zsize());
    mask_in_corr   = correct_NaN(mask);
    pv_map_in_corr = correct_NaN(pv_map);

    // Do correction on each slice of time series
    for(int i = 0; i < ndata_in; i++) {
      // Correct NaN and INF values of the 3D matrix of current TI (time domain)
      volume<float> corrected_data_ti = correct_NaN(data_in[i]);

      // Linear regression PV correction
      data_pvcorr[i] = correct_pv_lr(corrected_data_ti, mask_in_corr, pv_map_in_corr, kernel);
    }
  }

  // Function to correct PV using LR method
  volume<float> correct_pv_lr(const volume<float>& data_in, const volume<float>& mask, const volume<float>& pv_map, int kernel)
  {

    volume<float> submask;
    volume<float> data_roi;
    volume<float> pv_roi;
    RowVector pseudo_inv;
    Matrix pv_corr_result;

    // Variables to store the boundary index of submask (ROI)
    int x_0;
    int x_1;
    int y_0;
    int y_1;
    int z_0;
    int z_1;

    int count;
    float pv_ave = 0.0f;

    // Get x y z dimension
    int x = data_in.xsize();
    int y = data_in.ysize();
    int z = data_in.zsize();

    volume<float> corr_data(x, y, z); // result matrix

    // Linear regression to correct (smooth) the data
    for (int i = 0; i < x; i++) {
      for (int j = 0; j < y; j++) {
        for (int k = 0; k < z; k++) {
          // Only work with positive voxels
          if(mask.value(i, j, k) > 0) {

            // Determine ROI boundary index
            x_0 = max(i - kernel, 0);
            x_1 = min(i + kernel, x - 1);
            y_0 = max(j - kernel, 0);
            y_1 = min(j + kernel, y - 1);
            z_0 = max(k - kernel, 0);
            z_1 = min(k + kernel, z - 1);

            // create a submask here
            mask.setROIlimits(x_0, x_1, y_0, y_1, z_0, z_1);
            mask.activateROI();
            submask = mask.ROI();

            // calculate the sum of all elements in submask
            // proceed if sum is greater than 5 (arbitrary threshold)
            if(submask.sum() > 5) {
              /* Create an ROI (sub volume of data and PV map),
                then mask it with submask to create sub data and PV map */

              // Obtain ROI volume (must set limits and activate first)
              data_in.setROIlimits(x_0, x_1, y_0, y_1, z_0, z_1);
              pv_map.setROIlimits(x_0, x_1, y_0, y_1, z_0, z_1);
              data_in.activateROI();
              pv_map.activateROI();
              data_roi = data_in.ROI();
              pv_roi = pv_map.ROI();

              cout << "correctionhhh" << endl;
              getchar();
              
              
              ColumnVector data_roi_v_t = ColumnVector(data_roi.xsize() * data_roi.ysize() * data_roi.zsize());
              ColumnVector pv_roi_v_t = ColumnVector(pv_roi.xsize() * pv_roi.ysize() * pv_roi.zsize());

              count = 0;
              // Apply a mask on data_roi and pv_roi
              // Extract values from data_roi and pv_roi whose mask values are non-zero
              for(int a = 0; a < data_roi.xsize(); a++) {
                for(int b = 0; b < data_roi.ysize(); b++) {
                  for(int c = 0; c < data_roi.zsize(); c++) {
                    if(submask.value(a, b, c) > 0) {
                      data_roi_v_t.element(count) = data_roi.value(a, b, c);
                      pv_roi_v_t.element(count) = pv_roi.value(a, b, c);
                      count++;

                    }
                    else {
                      continue;
                    }
                  }
                }
              }
              
              ColumnVector data_roi_v = ColumnVector(count);
              ColumnVector pv_roi_v = ColumnVector(count);
              
              
              for(int a = 0; a < count; a++) {
                data_roi_v.element(a) = data_roi_v_t.element(a);
                pv_roi_v.element(a) = pv_roi_v_t.element(a);
              }
              
              // Unstable function
              //data_roi_v = apply_mask(data_roi, submask);
              //pv_roi_v = apply_mask(pv_roi, submask);

              
              // Now data_roi_v and pv_roi_v stores the non-zero elemnts within this submask

              // Deactivate ROI
              data_in.deactivateROI();
              pv_map.deactivateROI();

              // If pv_roi is all zeros, then the pseudo inversion matrix will be singular
              // This will cause run time error
              // So we assign the corrected result to zero in such cases
              if(pv_roi_v.IsZero()) {
                corr_data.value(i, j, k) = 0.0f;
              }
              else {
                // Compute pseudo inversion matrix of PV map
                // ((P^t * P) ^ -1) * (P^t)
                pseudo_inv = ( (pv_roi_v.t() * pv_roi_v).i() ) * (pv_roi_v.t());

                //cout << "  hh" << endl;

                // Get average PV value of the current kernel
                pv_ave = (float) pv_roi_v.Sum() / pv_roi_v.Nrows();

                // Calculate PV corrected data only if there is some PV compoment
                // If there is little PV small, make it zero
                if(pv_ave >= 0.01) {
                  pv_corr_result = pseudo_inv * data_roi_v;
                  corr_data.value(i, j, k) = pv_corr_result.element(0, 0);
                }
                else {
                  corr_data.value(i, j, k) = 0.0f;
                }
              }
            }

            else {
              // do nothing at the moment
            } // end submask

            // Discard current submask (ROI)
            mask.deactivateROI();

          }

          else{
            // do nothing at the moment
          } // end mask

        }
      }
    }

    return corr_data;

  } // End function correct_pv_lr

  // Function to correct NaN values
  volume<float> correct_NaN(const volume<float>& data_in) {

    // Clone the input data to output data
    volume<float> data_out = data_in;

    for(int i = 0; i < data_in.xsize(); i++) {
      for(int j = 0; j < data_in.ysize(); j++) {
        for(int k = 0; k < data_in.zsize(); k++) {
          // IEEE standard: comparison between NaN values is always false
          // i.e. NaN == NaN is false
          // In this case, we set it to zero
          if(data_in.value(i, j, k) != data_in.value(i, j, k)) {
            data_out.value(i, j, k) = 0.0f;
          }
          else {
            continue;
          }
        }
      }
    }

    return data_out;

  }





}
