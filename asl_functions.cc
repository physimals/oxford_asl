/*   asl_functions.cc various functions for the manipulation of ASL data

    Michael Chappell - FMRIB Image Analysis Group

    Moss Zhao - IBME Quantitative Biomedical Inference (QuBIc) Group

    Copyright (C) 2015 University of Oxford  */

/*   CCOPYRIGHT   */

#include "asl_functions.h"

namespace OXASL {

  void data2stdform(Matrix& datamtx, vector<Matrix>& asldata, int ntis, bool isblocked, bool ispairs, bool blockpairs) {
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
		  if (blockpairs) {
		    //we get all the tags for this repeat then all the controls (or vice versa)
		    thisti.Row(i) = datamtx.Row(2*(i-1)*ntis+ti);
		    thisti.Row(i+nrpts) = datamtx.Row(2*(i-1)*ntis+ntis+ti);
		    // NOTE that we keep it internally in blockpair format
		  }
		  else {
		    //tag controla pairs are adjacent volumes
		    thisti.Row(2*i-1) = datamtx.Row(2*(i-1)*ntis+2*ti-1);
		    thisti.Row(2*i) = datamtx.Row(2*(i-1)*ntis+2*ti);
		  }
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

  void separatepairs(vector<Matrix>& asldata, vector<Matrix>& asldataodd, vector<Matrix>& asldataeven, bool blockpairs) {
    int ntis = asldata.size();
    int nmeas = asldata[0].Nrows();
    int nrpts=nmeas/2; //if we are using this function then the data must contain pairs

    // just in case the vectors are not empty to start with
    asldataodd.clear();
    asldataeven.clear();

    int idx;
    
    for (int ti=0; ti<ntis; ti++) {
      
      Matrix oddmtx;
     
      Matrix evenmtx;
     

      if (blockpairs) {
	// data is in the form of a block of all tag (control) followed by a block of control (tag)
	oddmtx = asldata[ti].Rows(1,nrpts);
	evenmtx = asldata[ti].Rows(nrpts+1,nmeas);
	
      }
      else {
	// data is in the form of adjacent tag control pairs
	oddmtx = asldata[ti].Row(1);
	evenmtx = asldata[ti].Row(2);

	for (int r=2; r<=nrpts; r++) {
	  idx=(r-1)*2+1;
	  oddmtx &= asldata[ti].Row(idx);
	  evenmtx &= asldata[ti].Row(idx+1);
	}
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
	separatepairs(asldata,asldataodd,asldataeven,false);

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
      separatepairs(asldata,asldataodd,asldataeven,false);

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

  // Function to correct PV using LR method
  volume<float> correct_pv_lr(const volume<float>& data_in, const volume<float>& mask, const volume<float>& pv_map_gm, const volume<float>& pv_map_wm, int kernel)
  {

    volume<float> submask;
    volume<float> data_roi;
    volume<float> pv_roi;
    Matrix pseudo_inv;
    Matrix pv_corr_result;
    Matrix ha_result;

    int singular_matrix_flag = -1;

    // Variables to store the boundary index of submask (ROI)
    int x_0;
    int x_1;
    int y_0;
    int y_1;
    int z_0;
    int z_1;

    float pv_average = 0.0f;

    // Get x y z dimension
    int x = mask.xsize();
    int y = mask.ysize();
    int z = mask.zsize();

    volume<float> gm_corr_data(x, y, z); // result matrix GM
    //volume<float> wm_corr_data(x, y, z); // result matrix WM

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
            //mask.setROIlimits(x_0, x_1, y_0, y_1, z_0, z_1);
            //mask.activateROI();
            //submask = mask.ROI();

            // Define three column vectors to store data and PVE of the current regression kernel
            ColumnVector sub_mask = ColumnVector((x_1 - x_0 + 1) * (y_1 - y_0 + 1) * (z_1 - z_0 + 1));
            ColumnVector sub_data = ColumnVector((x_1 - x_0 + 1) * (y_1 - y_0 + 1) * (z_1 - z_0 + 1));
            Matrix sub_pve  = Matrix((x_1 - x_0 + 1) * (y_1 - y_0 + 1) * (z_1 - z_0 + 1), 2);
            
            int sub_mask_count = 0;
            int non_zero_count = 0;
            float submask_sum = 0.0f; // value to store the sum of current mask kernel
            for(int p = 0; p < z_1 - z_0 + 1; p++) {
              for(int n = 0; n < y_1 - y_0 + 1; n++) {
                for(int m = 0; m < x_1 - x_0 + 1; m++) {
                  //sub_mask.element(sub_mask_count) = mask.value(x_0 + m, y_0 + n, z_0 + p);
                  //sub_data.element(sub_mask_count) = data_in.value(x_0 + m, y_0 + n, z_0 + p);
                  //sub_pve.element(sub_mask_count) = pv_map.value(x_0 + m, y_0 + n, z_0 + p);
                  sub_mask.element(sub_mask_count) = mask.value(x_0 + m, y_0 + n, z_0 + p);
                  sub_data.element(sub_mask_count) = data_in.value(x_0 + m, y_0 + n, z_0 + p);
                  // In the Sub PVE matrix, first column is GM
                  sub_pve.element(sub_mask_count, 0) = pv_map_gm.value(x_0 + m, y_0 + n, z_0 + p);
                  // In the Sub PVE matrix, second column is WM
                  sub_pve.element(sub_mask_count, 1) = pv_map_wm.value(x_0 + m, y_0 + n, z_0 + p);
                  submask_sum = submask_sum + mask.value(x_0 + m, y_0 + n, z_0 + p);
                  if(mask.value(x_0 + m, y_0 + n, z_0 + p) > 0) {
                    non_zero_count++;
                  }
                  sub_mask_count++;
                }
              }
            }

            // calculate the sum of all elements in submask
            // proceed if sum is greater than 5 (arbitrary threshold)
            if(submask_sum > 5) {
              // Apply submask to the data and PVE of the current kernel
              ColumnVector data_roi_v = ColumnVector(non_zero_count);
              Matrix pv_roi_v = Matrix(non_zero_count, 2);
              //RowVector pv_roi_r = RowVector(non_zero_count);
              
              int non_zero_index = 0;
              // Extract all non-zero elements
              //cout << sub_mask_count << endl;
              for(int a = 0; a < sub_mask_count; a++) {
                if(sub_mask.element(a) > 0) {
                  data_roi_v.element(non_zero_index) = sub_data.element(a);
                  pv_roi_v.element(non_zero_index, 0) = sub_pve.element(a, 0);
                  pv_roi_v.element(non_zero_index, 1) = sub_pve.element(a, 1);
                  non_zero_index++;
                }
                else {
                  continue;
                }
              }

              // If pv_roi is all zeros, then the pseudo inversion matrix will be singular
              // This will cause run time error
              // So we assign the corrected result to zero in such cases
              float det = ((pv_roi_v.t()) * pv_roi_v).Determinant();
              if( (det <= 0.00001) && (det >= 0) ) {
                gm_corr_data.value(i, j, k) = 0.0f;
                singular_matrix_flag = 0;
                //cout << i << ", " << j << ", " << k << endl;
                //cout << "singular" << endl;
                //getchar();
              }
              else {
                // Compute pseudo inversion matrix of PV map
                // ((P^t * P) ^ -1) * (P^t)
                //float haha = ((pv_roi_v.t()) * pv_roi_v).Determinant();
                //cout << ((pv_roi_v.t()) * pv_roi_v).Determinant() << endl;
                pseudo_inv = ( ( (pv_roi_v.t()) * pv_roi_v).i() ) * (pv_roi_v.t());
                //cout << i << ", " << j << ", " << k << endl;
                // Get average PV value of the current kernel
                pv_average = (float) pv_roi_v.Sum() / pv_roi_v.Nrows();

                // Calculate PV corrected data only if there is some PV compoment
                // If there is little PV small, make it zero
                if(pv_average >= 0.01) {
                  pv_corr_result = pseudo_inv * data_roi_v;
                  gm_corr_data.value(i, j, k) = pv_corr_result.element(0, 0); // output GM only
                  //corr_data.value(i, j, k) = pv_corr_result.element(1, 0);
                }
                else {
                  gm_corr_data.value(i, j, k) = 0.0f;
                }

              }
            }

            else {
            } // end submask

          }

          else{
            // do nothing at the moment
          } // end mask

        }
      }
    }

    if(singular_matrix_flag == 0) {
      cout << "Caution: singular matrix found in PV Correction. This usually happens to data from Siemens. No action required." << endl;
    }

    return gm_corr_data;

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


  // function to perform partial volume correction by linear regression
  void pvcorr_LR(const volume4D<float>& data_in, int ndata_in, const volume<float>& mask, const volume<float>& pv_map_gm, const volume<float>& pv_map_wm, int kernel, volume4D<float>& data_pvcorr) {

    // Version control
    cout << "PV correction by linear regression. version 1.0.4 (beta). Last compiled on 20160203" << endl;
    
    // Clone input data to pv corrected data
    data_pvcorr = data_in;

    // Correct NaN and INF numbers of input mask and pvmap
    volume<float> mask_in_corr(mask.xsize(), mask.ysize(), mask.zsize());
    volume<float> pv_map_gm_in_corr(pv_map_gm.xsize(), pv_map_gm.ysize(), pv_map_gm.zsize());
    volume<float> pv_map_wm_in_corr(pv_map_gm.xsize(), pv_map_gm.ysize(), pv_map_gm.zsize());
    mask_in_corr   = correct_NaN(mask);
    pv_map_gm_in_corr = correct_NaN(pv_map_gm);
    pv_map_wm_in_corr = correct_NaN(pv_map_wm);

    // Do correction on each slice of time series
    for(int i = 0; i < ndata_in; i++) {
      // Correct NaN and INF values of the 3D matrix of current TI (time domain)
      volume<float> corrected_data_ti = correct_NaN(data_in[i]);

      // Linear regression PV correction
      data_pvcorr[i] = correct_pv_lr(corrected_data_ti, mask_in_corr, pv_map_gm_in_corr, pv_map_wm_in_corr, kernel);
    }
  }


  // function to convert PAR REC file to NIfTI format
  void convert_par_rec_to_nifti(const string file_par_name, const string file_rec_name, volume<float>& mask_nifti, volume4D<float>& data_nifti) {

    // variables for 4D matrix position index
    // In total, there should be 49 variables
    int variables_count = 49;

    int slice_number;
    int echo_number;
    int dynamic_scan_number;
    int cardiac_phase_number;
    int image_type_mr;
    int scanning_sequence;
    int index_in_REC_file; // in images
    int image_pixel_size; // in bits
    int scan_percentage;
    int recon_resolution_x;
    int recon_resolution_y;
    float rescale_intercept;
    float rescale_slope;
    float scale_slope;
    int window_center;
    int window_width;
    float image_angulation_ap; // in degrees
    float image_angulation_fh; // in degrees
    float image_angulation_rl; // in degrees
    float image_offcentre_ap; // in mm
    float image_offcentre_fh; // in mm
    float image_offcentre_rl; //in mm
    float slice_thickness; // in mm
    float slice_gap; // in mm
    int image_display_orientation;
    int slice_orientation; // TRA/SAG/COR
    int fmri_status_indication;
    int image_type_ed_es; // end diast/end syst
    float pixel_spacing_x; // in mm
    float pixel_spacing_y; // in mm
    float echo_time;
    float dyn_scan_begin_time;
    float trigger_time;
    float diffusion_b_factor;
    int number_of_averages;
    float image_flip_angle; // in degrees
    int cardiac_frequency; // bm;
    int minimum_RR_interval; // in ms
    int maximum_RR_interval; // in ms
    int turbo_factor; // 0=no turbo
    float inversion_delay; // in ms
    int diffusion_b_value_number; // imagekey!
    int gradient_orientation_number; // imagekey!
    string constrast_type;
    string diffusion_anisotropy_type;
    float diffusion_ap;
    float diffusion_fh;
    float diffusion_rl;
    int label_type; // imagekey, for ASL. 1=control; 2=label

    // start file handling
    ifstream par_file(file_par_name.c_str());
    ifstream rec_file(file_rec_name.c_str(), ifstream::binary);
    int read_counter = 0;

    int x_y_dimension = data_nifti.xsize();
    short buffer[x_y_dimension * x_y_dimension];
    float voxel_floating_point_value[x_y_dimension * x_y_dimension];

    // define a temporary matrix [x][y][z][dynamics (phases)][cardiac phases (TIs)][ASL (value: 1 control; 2 control)]
    //int x_dim = 64, y_dim = 64, z_dim = 15, dynamics = 7, shift = 2, repeats = 1, tis = 11, asl = 2; // ASL file
    int x_dim = 64, y_dim = 64, z_dim = 15, dynamics = 7, shift = 2, repeats = 2, tis = 11, asl = 2; // ASL file
    //int x_dim = 288, y_dim = 288, z_dim = 245, dynamics = 1, shift = 1, repeats = 1, tis = 1, asl = 1; // Structure file

    // Display dimension of the output file
    cout << "x: " << x_dim << " y: " << y_dim << " z: " << z_dim << " t: " << dynamics * shift * repeats * tis * asl << endl;

    static float temp_data[64][64][15][7 * 2 * 2][11][2];
    //static float temp_data[288][288][245][1][1][1];
    //static float temp_data_2[64][64][15][7][11 * 2][2]; // ASL file array
    //int structural_image_flag = 1;
    int structural_image_flag = 0;

    string current_line;
    int line_count = 0;
    int data_start_line_number = 101;
    int data_end_line_number = data_start_line_number + z_dim * dynamics * shift * repeats * tis * asl - 1; // ASL file single repeats
    
    /*
    // Check big or small endian
    int num = 1;
    if(*(char *)&num == 1) {
      printf("\nLittle-Endian\n");
    }
    else {
        printf("Big-Endian\n");
    }
    */
    
    if(par_file.is_open()) {

      while(getline(par_file, current_line)) {
        
        line_count++;

        // Empty line
        if(current_line == "") {
          continue;
        }

        // Comment line (start with #)
        if(current_line[0] == '#') {
          continue;
        }

        // Metadata line (start with .)
        if(current_line[0] == '.') {
          continue;
        }

        // position to start copying data
        if(line_count >= data_start_line_number && line_count <= data_end_line_number) {
        //if(line_count >= 101 && line_count <= 4720) {
          float current_line_data[variables_count];

          stringstream current_line_in_stream(current_line);

          // split current line's data into elements in the array
          for(int i = 0; i < variables_count; i++) {
            current_line_in_stream >> current_line_data[i];
            //cout << current_line_data[i] << endl;
            //getchar();
          }

          // Variable to deal with converting float to strings
          ostringstream str_buff;

          // Now current_line_data saves each number in the current line in floating point format
          slice_number                = (int) current_line_data[0];
          echo_number                 = (int) current_line_data[1];
          dynamic_scan_number         = (int) current_line_data[2];
          cardiac_phase_number        = (int) current_line_data[3];
          image_type_mr               = (int) current_line_data[4];
          scanning_sequence           = (int) current_line_data[5];
          index_in_REC_file           = (int) current_line_data[6]; // in images
          image_pixel_size            = (int) current_line_data[7]; // in bits
          scan_percentage             = (int) current_line_data[8];
          recon_resolution_x          = (int) current_line_data[9];
          recon_resolution_y          = (int) current_line_data[10];
          rescale_intercept           = current_line_data[11];
          rescale_slope               = current_line_data[12];
          scale_slope                 = current_line_data[13];
          window_center               = (int) current_line_data[14];
          window_width                = (int) current_line_data[15];
          image_angulation_ap         = current_line_data[16]; // in degrees
          image_angulation_fh         = current_line_data[17]; // in degrees
          image_angulation_rl         = current_line_data[18]; // in degrees
          image_offcentre_ap          = current_line_data[19]; // in mm
          image_offcentre_fh          = current_line_data[20]; // in mm
          image_offcentre_rl          = current_line_data[21]; //in mm
          slice_thickness             = current_line_data[22]; // in mm
          slice_gap                   = current_line_data[23]; // in mm
          image_display_orientation   = (int) current_line_data[24];
          slice_orientation           = (int) current_line_data[25]; // TRA/SAG/COR
          fmri_status_indication      = (int) current_line_data[26];
          image_type_ed_es            = (int) current_line_data[27]; // end diast/end syst
          pixel_spacing_x             = current_line_data[28]; // in mm
          pixel_spacing_y             = current_line_data[29]; // in mm
          echo_time                   = current_line_data[30];
          dyn_scan_begin_time         = current_line_data[31];
          trigger_time                = current_line_data[32];
          diffusion_b_factor          = current_line_data[33];
          number_of_averages          = (int) current_line_data[34];
          image_flip_angle            = current_line_data[35]; // in degrees
          cardiac_frequency           = (int) current_line_data[36]; // bm;
          minimum_RR_interval         = (int) current_line_data[37]; // in ms
          maximum_RR_interval         = (int) current_line_data[38]; // in ms
          turbo_factor                = (int) current_line_data[39]; // 0=no turbo
          inversion_delay             = current_line_data[40]; // in ms
          diffusion_b_value_number    = (int) current_line_data[41]; // imagekey!
          gradient_orientation_number = (int) current_line_data[42]; // imagekey!
          str_buff << current_line_data[43];
          constrast_type              = str_buff.str();
          str_buff << current_line_data[44];
          diffusion_anisotropy_type   = str_buff.str();
          //constrast_type              = to_string(current_line_data[43]);
          //diffusion_anisotropy_type   = to_string(current_line_data[44]);
          diffusion_ap                = current_line_data[45];
          diffusion_fh                = current_line_data[46];
          diffusion_rl                = current_line_data[47];
          label_type                  = (int) current_line_data[48]; // imagekey, for ASL. 1=control; 2=label

          // Now read 16bit integer from REC file
          // We need to read 64x64 numbers (for each slice)
          if(rec_file.read((char*)&buffer, sizeof(buffer))) {
            read_counter++;
          }

          // calculate the corresponding floating point value
          for(int k = 0; k < (sizeof(buffer) / sizeof(buffer[0])); k++) {
            voxel_floating_point_value[k] = convert_pixel_value_to_floating_point((int)buffer[k], rescale_slope, rescale_intercept, scale_slope);
          }

          // now we need to assign this to a position

          // 22(TIs) * 7(Phases) * 2(Tag-control)
          // Hierarchy of the t dimension of nifti file (from top to bottom):
          // 22 chuncks, each containing TIs
          // In each chuck, there are 7 phases
          // In each phase, there are the tag-control pair

          // useful variables here: 
          // slice_number: which z index (1-15)
          // dynamic_scan_number: which phase (dynamic 1: 1-7; dynamic 2: 8-14)
          // cardiac_phase_number: which TI (1-11)
          // label_type: tag or control (1 control, 2 tag)

          //int z_index = slice_number - 1;
          //int t_index = dynamic_scan_number * cardiac_phase_number * label_type - 1; // index start from zero, so we need to subtract one

          // define a temporary matrix [x][y][z][dynamics (phases)][cardiac phases (TIs)][ASL (value: 1 control; 2 control)]
          // int x_dim = 64, y_dim = 64, z_dim = 15, dynamics = 7, shift = 2, tis = 11, asl = 2;
          // float temp_data[x_dim][y_dim][z_dim][dynamics * shift][tis][asl];
          // float temp_data_2[x_dim][y_dim][z_dim][dynamics][tis * shift][asl];

          int n = 0;
          for(int x = 0; x < data_nifti.xsize(); x++) {
            for(int y = 0; y < data_nifti.ysize(); y++) {

              // case for structural image
              if(structural_image_flag == 1) {
                data_nifti.value(x, y, slice_number - 1, 0) = voxel_floating_point_value[n];
              }
              temp_data[x][y][slice_number - 1][dynamic_scan_number - 1][cardiac_phase_number - 1][label_type - 1] = voxel_floating_point_value[n];
              ++n;
            }
          }

        }


      }

      par_file.close();
      rec_file.close();
    }
    
    // Put the data from array to output 4D matrix
    int t = 0;
    for(int x = 0; x < x_dim; x++) {
      for(int y = 0; y < y_dim; y++) {
        for(int z = 0; z < z_dim; z++) {

          for(int a = 0; a < dynamics * shift * repeats; a++) {
            for(int b = 0; b < tis; b++) {
              for(int c = 0; c < asl; c++) {
                data_nifti.value(x, y, z, t) = temp_data[x][y][z][a][b][c];
                t++;
              }
            }
          }

          t = 0;
        }
      }
    }
    
    // End of function

  }


  // function to make default mask
  void create_default_mask(volume<float>& mask_nifti) {
    // make a mask
    float mask_intensity = 10.0f; // arbitrary value
    
    for(int x = 0; x < mask_nifti.xsize(); x++) {
      for(int y = 0; y < mask_nifti.ysize(); y++) {
        for(int z = 0; z < mask_nifti.zsize(); z++) {
          mask_nifti.value(x, y, z) = mask_intensity;
        }
      }
    }

  }

  // function to make default data nifti file
  void create_default_data_nifti(volume4D<float>& data_nifti) {

    // make a default 4D data file
    float voxel_intensity = -999998.0f; // arbitrary value
    
    for(int x = 0; x < data_nifti.xsize(); x++) {
      for(int y = 0; y < data_nifti.ysize(); y++) {
        for(int z = 0; z < data_nifti.zsize(); z++) {
          for(int t = 0; t < data_nifti.tsize(); t++) {
            data_nifti.value(x, y, z, t) = voxel_intensity;
          }
        }
      }
    }

  }

  float convert_pixel_value_to_floating_point(int pixel_value, float rescale_slope, float rescale_intercept, float scale_slope) {

    float displayed_value_on_console = (float) pixel_value * rescale_slope + rescale_intercept;

    float floating_point_value = displayed_value_on_console / (rescale_slope * scale_slope);

    return floating_point_value;

  }

  Matrix extrapolate_avg(Matrix data_in, Matrix mask_in, int neighbour_size) {

    Matrix data_extrapolated = data_in;

    int x = data_in.Ncols();
    int y = data_in.Nrows();

    int x_index = 0;
    int y_index = 0;
    int dx = 0;
    int dy = -1;

    int x_boundary = x - 1;
    int y_boundary = y - 1;
    int x_offset = x / 2;
    int y_offset = y / 2;

    int t = max(x_boundary, y_boundary);
    int max_i = t * t;

    int count = 1;

    for (int i = 0; i < max_i; i++) {

      // Position found
      if ( (-x_boundary / 2 <= x_index) && (x_index <= x_boundary / 2) && (-1 * y_boundary / 2 <= y_index) && (y_index <= y_boundary / 2)) {
        //cout << x_index << ", " << y_index << endl;

        // Do extrapolation
        int x_index_on_matrix = x_index + x_offset;
        int y_index_on_matrix = y_index + y_offset;

        // Only work on eroded voxels
        if (mask_in.element(x_index_on_matrix, y_index_on_matrix) != 0 && data_in.element(x_index_on_matrix, y_index_on_matrix) == 0 ){
          // Create a square matrix of size neighbourhood and centered at the current postion
          int off_set = floor(neighbour_size / 2);

          int column_begin = x_index_on_matrix - off_set;
          int column_end = x_index_on_matrix + off_set;

          int row_begin = y_index_on_matrix - off_set;
          int row_end = y_index_on_matrix + off_set;

          // If it is out of boundary then continue
          if (column_begin <= 0 || column_end <= 0 || row_begin <= 0 || row_end <= 0) {
            continue;
          }

          float sum = 0;
          int non_zero_count = 0;

          for (int m = column_begin; m <= column_end; m++) {
            for (int n = row_begin; n <= row_end; n++) {
              if (data_in.element(m, n) != 0) {
                sum = sum + data_in.element(m, n);
                non_zero_count++;
              }
            }
          }

          if(non_zero_count > 0) {
            data_extrapolated.element(x_index_on_matrix, y_index_on_matrix) = sum / non_zero_count;
          }

        }

      }

      if( (x_index == y_index) || ((x_index < 0) && (x_index == (-1) * y_index)) || ((x_index > 0) && (x_index == 1 - y_index))  ) {
        t = dx;
        dx = -1 * dy;
        dy = t;
      }

      x_index = x_index + dx;
      y_index = y_index + dy;
    }

    return data_extrapolated;

  }

  // Function to extrapolate voxels
  void extrapolate(const volume4D<float>& data, int ndata_in, const volume<float>& mask, int neighbour_size, volume4D<float>& data_extrapolated) {

    // Version control
    cout << "Extrapolation. version 1.0.1 (beta). Last compiled on 20161029" << endl;
    
    // Clone input data to pv corrected data
    data_extrapolated = data;

    // Correct NaN and INF numbers of input mask and pvmap
    volume<float> mask_in_corr(mask.xsize(), mask.ysize(), mask.zsize());
    //volume<float> pv_map_gm_in_corr(pv_map_gm.xsize(), pv_map_gm.ysize(), pv_map_gm.zsize());
    //volume<float> pv_map_wm_in_corr(pv_map_gm.xsize(), pv_map_gm.ysize(), pv_map_gm.zsize());
    mask_in_corr = correct_NaN(mask);
    //pv_map_gm_in_corr = correct_NaN(pv_map_gm);
    //pv_map_wm_in_corr = correct_NaN(pv_map_wm);

    // Do correction on each slice of time series
    for(int i = 0; i < ndata_in; i++) {
      // Correct NaN and INF values of the 3D matrix of current TI (time domain)
      volume<float> nan_corrected_data_ti = correct_NaN(data[i]);

      // Define a temporary matrix to save current extrapolated results
      volume<float> extrapolated_data_3D(mask.xsize(), mask.ysize(), mask.zsize());

      // Get x y z dimension
      int x = nan_corrected_data_ti.xsize();
      int y = nan_corrected_data_ti.ysize();
      int z = nan_corrected_data_ti.zsize();

      // for each slice, perform extrapolation
      for(int j = 0; j < z; j++) {
        Matrix data_slice_non_extrapolated = Matrix(x, y);
        Matrix data_slice_extrapolated = Matrix(x, y);
        Matrix data_mask = Matrix(x, y);

        // Copy the current slice to non-extrapolated matrix
        for(int m = 0; m < x; m++) {
          for(int n = 0; n < y; n++) {
            data_slice_non_extrapolated.element(m, n) = nan_corrected_data_ti.value(m, n, j);
            data_mask.element(m, n) = mask_in_corr.value(m, n, j);
            // Default value for the extrapolated matrix is the same with the input file
            data_slice_extrapolated.element(m, n) = nan_corrected_data_ti.value(m, n, j);
          }
        }

        data_slice_extrapolated = extrapolate_avg(data_slice_non_extrapolated, data_mask, neighbour_size);

        // Assign the extrapolated matrix to the 3D volume
        for(int m = 0; m < x; m++) {
          for(int n = 0; n < y; n++) {
           extrapolated_data_3D.value(m, n, j) = data_slice_extrapolated.element(m, n);
          }
        }


      }

      // Assign result to the 4D volume
      data_extrapolated[i] = extrapolated_data_3D;
    }
  }

}
