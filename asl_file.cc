/*   asl_file.cc file manipulator for multi-TI ASL data

    Michael Chappell - FMRIB Image Analysis Group

    Moss Zhao - IBME Quantitative Biomedical Inference (QuBIc) Group

    Copyright (C) 2015 University of Oxford  */

/*   CCOPYRIGHT   */

#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"
#include "newmatap.h"
#include "newmatio.h"
#include "stdlib.h"
#include "utils/tracer_plus.h"
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>

#include "asl_functions.h"
#include "readoptions.h"

using namespace Utilities;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace MISCMATHS;

using namespace OXASL;

int main(int argc, char *argv[])
{
    try
    {
        //parse command line (puts these into the log file)
        ReadOptions &opts = ReadOptions::getInstance();
        opts.parse_command_line(argc, argv);

        //deal with input data type options
        bool isblocked = false;  //indicates if data is in blocks of repeats (of all TIs) rather than TIs
        bool ispairs = false;    //indicates if data contains adjacent pairs of measurments
        bool isdiff = false;     //indicates if we have differenced data
        bool tagfirst = true;    //indicates that tag comes first in tag-control pairs
        bool blockpairs = false; //indicates that the tag-control pairs are grouped separelty within the blocks, not in adacent volumes of tag then control
        //ispairs=opts.ispairs.value();

        // block format: isblocked indicates if data are in *blocks of repeats*
        string ibf = opts.inblockform.value();
        if (ibf.compare("rpt") == 0)
            isblocked = true;
        else if (ibf.compare("tis") == 0)
            isblocked = false;
        else
            throw Exception("Unrecognised input block format");

        string iaf = opts.inaslform.value();
        if (iaf.compare("diff") == 0)
            isdiff = true; // differenced data
        else if (iaf.compare("tc") == 0)
            isdiff = false; //tag control pairs, tag first
        else if (iaf.compare("ct") == 0)
        {
            isdiff = false;
            tagfirst = false;
        } //tag control pairs, control first
        else if (iaf.compare("tcb") == 0)
        {
            isdiff = false;
            blockpairs = true;
        } // tag control pairs, tag first, all the tags and controls are sepratedly grouped together (not adjacent to each other)
        else if (iaf.compare("ctb") == 0)
        {
            isdiff = false;
            blockpairs = true;
            tagfirst = false;
        } // ditto, control first
        else
            throw Exception("Unrecognised input asl form");

        ispairs = !isdiff; //convienient to have ispairs as a bool

        bool outblocked;
        string obf = opts.outblockform.value();
        if (obf.compare("rpt") == 0)
            outblocked = true;
        else if (obf.compare("tis") == 0)
            outblocked = false;
        else if (obf.compare("notset") == 0)
            outblocked = isblocked; //default is to set output same as input
        else
            throw Exception("Unrecognised output block format");

        /*
    bool outdiff;
    string oaf=opts.outaslform.value();
    if     (oaf.compare("diff")==0) outdiff=true;
    else if(oaf.compare("tc")==0) outdiff=false;
    else throw Exception("Unrecognised output asl form");
    */
        bool outpairs = ispairs; // outpairs indicates wehter the data we are processing for output is in the form of pairs - by deafult if input is in pairs then output the pairs

        //load data
        volume4D<float> data;
        read_volume4D(data, opts.datafile.value());

        // Partail volume correction variables
        volume<float> pv_gm_map;
        volume<float> pv_wm_map;
        int kernel = 0;
        vector<Matrix> output_data_gm;
        vector<Matrix> output_data_wm;
        // Read partial volume map and kernel size
        if (opts.pv_gm_file.set() && opts.pv_wm_file.set() && opts.kernel.set())
        {
            read_volume(pv_gm_map, opts.pv_gm_file.value());
            read_volume(pv_wm_map, opts.pv_wm_file.value());
            kernel = opts.kernel.value();
        }

        // Extrapolation variables
        // To save output data
        vector<Matrix> output_data_extrapolated;
        // Extrapolation options
        int neighbour_size = 0;
        if (opts.extrapolate_option.value())
        {
            // Get the neighbourhood size
            neighbour_size = opts.neighbour.value();
        }

        // load mask
        // if a mask is not supplied then default to processing whole volume
        // Create a mask from the data volume to ensure metadata is correct
        volume<float> mask = data[0] * 0;
        mask = 1;
        if (opts.maskfile.set())
            read_volume(mask, opts.maskfile.value());

        Matrix datamtx;
        datamtx = data.matrix(mask);
        int nvox = datamtx.Ncols();

        cout << "Number of voxels is:" << nvox << endl;

        //establish number of repeats

        // deal with the repeats string if specified
        int ntis = opts.ntis.value();
        int ndata = 0;
        vector<int> nrpts;

        string rpts = opts.rpts.value();
        if (rpts.compare("NULL") == 0)
        {
            int nmeas = data.tsize() / ntis; //number of measurements at each TI
            if (ispairs)
                nrpts.push_back(nmeas / 2); // number of repeats, note that is data is pairs then we have half as many true repeats
            else
                nrpts.push_back(nmeas);
            cout << "Number of repeats in data is:" << nrpts[0] << endl;

            // calculate the expected size of the data
            if (ispairs)
                ndata = ntis * nrpts[0] * 2;
            else
                ndata = nmeas * ntis;
        }
        else
        {
            //we have a comma separated list of repeats
            if (isblocked)
                throw Exception("Cannot have different number of repeats at each TI when --ibf=rpt, do not specify --rpts this will be calcuateld from the data");
            vector<string> result;
            stringstream ss(rpts);
            while (ss.good())
            {
                string substr;
                getline(ss, substr, ',');
                stringstream ss(substr);
                int thisnrpts;
                ss >> thisnrpts;
                nrpts.push_back(thisnrpts);

                cout << "Number of repeats at TI" << nrpts.size() << " is " << thisnrpts << endl;

                //cumulative calculation of the expected size of the data
                if (ispairs)
                    ndata += thisnrpts * 2;
                else
                    ndata += thisnrpts;
            }
            if (int(nrpts.size()) < ntis)
                throw Exception("Insufficient entries in list fo repeats (--rpts) compared to number of TIs");
        }

        if (ndata < data.tsize())
        {
            //some leftover samples in the data, produce a warning
            cout << "Warning: spare measurements found at end of data!" << endl
                 << "    Number of measurements discarded from end is: " << data.tsize() - ndata << endl;

            // discard
            datamtx = datamtx.Rows(1, ndata);
        }

        //get data into 'standard' format
        vector<Matrix> asldata;
        data2stdform(datamtx, asldata, ntis, nrpts, isblocked, ispairs, blockpairs);

        //deal with the splitting of pairs
        vector<Matrix> asldataodd;
        vector<Matrix> asldataeven;
        // int idx;

        if (opts.splitpairs.value() && opts.tcdiff.value())
        {
            // doesn't make sense to try and do both splitpairs and tcdifference
            throw Exception("Cannot both split the pairs and difference them");
        }

        if (opts.splitpairs.value())
        {
            // need to split the data here if split pairs has been requested
            separatepairs(asldata, asldataodd, asldataeven);
        }

        //tag control difference
        if (opts.tcdiff.value())
        {
            separatepairs(asldata, asldataodd, asldataeven); //split pairs ready for differencing
            //overwrite asldata with differenced data
            for (int ti = 0; ti < ntis; ti++)
            {
                asldata[ti] = asldataeven[ti] - asldataodd[ti];
                if (!tagfirst)
                    asldata[ti] *= -1.0; //if control image is first then the sign will be wrong here
            }
            outpairs = false;
        }

        // surround tag-control difference
        if (opts.surrtcdiff.value())
        {
            separatepairs(asldata, asldataodd, asldataeven); //split pairs ready for differencing
            //overwrite asldata with differenced data
            for (int ti = 0; ti < ntis; ti++)
            {
                //Matrix tempindata;
                //tempindata=asldata[ti];
                Matrix tempdata(asldata[ti].Nrows() - 1, asldata[ti].Ncols()); //surr diff data will have one fewer measurement than non-differenced data
                for (int aq = 1; aq < tempdata.Nrows(); aq = aq + 2)
                {
                    tempdata.Row(aq) = asldataeven[ti].Row(aq) - asldataodd[ti].Row(aq);
                    tempdata.Row(aq + 1) = asldataeven[ti].Row(aq) - asldataodd[ti].Row(aq + 1);
                }
                asldata[ti] = tempdata;
                if (!tagfirst)
                    asldata[ti] *= -1.0; //if control image is first then the sign will be wrong here
            }
            outpairs = false;
        }

        // Here we throw some exceptions that we currently don't support
        // PVC and Split ASL file
        if ((opts.pv_gm_file.set() || opts.pv_wm_file.set()) && (opts.splitpairs.value()))
        {
            throw Exception("We don't support partial volume correction and split files done together!");
        }
        // PVC and Extrapolation
        if ((opts.pv_gm_file.set() || opts.pv_wm_file.set()) && (opts.extrapolate_option.value()))
        {
            throw Exception("We don't support partial volume correction and extrapolation done together!");
        }
        if ((opts.splitpairs.value()) && (opts.extrapolate_option.value()))
        {
            throw Exception("We don't support split files and extrapolation done together!");
        }

        // Partial Volume Correction Options
        // Partial volume correction on each TI
        if (opts.pv_gm_file.set() && opts.pv_wm_file.set())
        {
            // Check mask file is specified
            if ((opts.maskfile.set()) && (opts.kernel.set()) && (opts.out.set()))
            {
                cout << "Start partial volume correction" << endl;

                cout << "Dealing with GM PV Correction" << endl;
                pvcorr_LR(asldata, ndata, mask, pv_gm_map, pv_wm_map, kernel, output_data_gm, outblocked, outpairs, nrpts, isblocked, ispairs, blockpairs);

                cout << "Dealing with WM PV Correction" << endl;
                pvcorr_LR(asldata, ndata, mask, pv_wm_map, pv_gm_map, kernel, output_data_wm, outblocked, outpairs, nrpts, isblocked, ispairs, blockpairs);

                cout << "Partial volume correction done!" << endl;
            }
            else if (!opts.maskfile.set())
            {
                throw Exception("Missing mask file. --mask=<mask file>");
            }
            else if (!opts.kernel.set())
            {
                throw Exception("Missing kernel size. --kernel=<3 to 9 integer>");
            }
            else if (!opts.out.set())
            {
                throw Exception("Missing output file. --out=<output file name>");
            }
            else
            {
                throw Exception("Halt!");
            }
        }

        // Extrapolation options
        if ((opts.extrapolate_option.value()) && (opts.out.set()))
        {
            // Check mask and input file
            if (opts.maskfile.set())
            {
                cout << "Start extrapolation!" << endl;

                // Perform extrapolation, result in output_data_extrapolated
                extrapolate(asldata, ndata, mask, neighbour_size, output_data_extrapolated, outblocked, outpairs, nrpts, isblocked, ispairs, blockpairs);
            }

            else if (!opts.extrapolate_option.value())
            {
                throw Exception("Missing mask file. --mask=<mask file>");
            }
            else if (!opts.out.set())
            {
                throw Exception("Missing output file. --out=<out file>");
            }
            else
            {
                throw Exception("Halt!");
            }
        }

        vector<Matrix> asldataout; //the data we are about to use for output purposes
        int nout = 1;              //number of output cycles we need to go through
        string fsub = "";

        if (opts.splitpairs.value())
        {
            nout = 2;
            outpairs = false;
        }

        if (opts.pv_gm_file.set() && opts.pv_wm_file.set())
        {
            nout = 2;
        }

        // OUTPUT: Main output section - anything that is compatible with the splitting of pairs (splitpairs) can go in here
        for (int o = 1; o <= nout; o++)
        {
            if (!opts.splitpairs.value())
                asldataout = asldata;
            else if (o == 1)
            {
                asldataout = asldataodd;
                fsub = "_odd";
                cout << "Dealing with odd members of pairs" << endl;
            }
            else
            {
                asldataout = asldataeven;
                fsub = "_even";
                cout << "Dealing with even members of pairs" << endl;
            }

            // Partial volume correction on each TI
            if (opts.pv_gm_file.set() && opts.pv_wm_file.set())
            {
                /*
        // Check mask file is specified
        if( (opts.maskfile.set()) && (opts.kernel.set()) && (opts.out.set()) )  {

          cout << "Start partial volume correction" << endl;


          // Convert asldataout to volume4D<float>
          Matrix aslmatrix_non_pvcorr;
          volume4D<float> asldata_non_pvcorr;
          stdform2data(asldataout, aslmatrix_non_pvcorr, outblocked, outpairs);
          asldata_non_pvcorr.setmatrix(aslmatrix_non_pvcorr, mask);

          if(o == 1) {
            cout << "Dealing with GM PV Correction" << endl;
            fsub = "_gm";
            pvcorr_LR(asldata_non_pvcorr, ndata, mask, pv_gm_map, pv_wm_map, kernel, data_pvcorr);
          }

          if(o == 2) {
            cout << "Dealing with WM PV Correction" << endl;
            fsub = "_wm";
            pvcorr_LR(asldata_non_pvcorr, ndata, mask, pv_wm_map, pv_gm_map, kernel, data_pvcorr);
          }
          // function to perform partial volume correction by linear regression
          //pvcorr_LR(asldata_non_pvcorr, ndata, mask, pv_gm_map, pv_wm_map, kernel, data_pvcorr);

          //covert data_pvcorr to vector<Matrix> aka stdform 
          Matrix data_pvcorr_mtx;
          vector<Matrix> asldataout_pvcorr;
          data_pvcorr_mtx = data_pvcorr.matrix(mask);
          data2stdform(data_pvcorr_mtx, asldataout_pvcorr, ndata, nrpts, isblocked, ispairs,blockpairs);
          asldataout = asldataout_pvcorr;

          //save_volume4D(data_pvcorr, pvout_file_name);
          cout << "Partial volume correction done!" << endl;
        }
        else if(!opts.maskfile.set()) {
          throw Exception("Missing mask file. --mask=<mask file>");
        }
        else if(!opts.kernel.set()) {
          throw Exception("Missing kernel size. --kernel=<3 to 9 integer>");
        }
        else if(!opts.out.set()) {
          throw Exception("Missing output file. --out=<output file name>");
        }
        else {
          throw Exception("Halt!");
        }
        */

                // GM case
                if (o == 1)
                {
                    //cout << "Dealing with GM PV Correction" << endl;
                    fsub = "_gm";
                    asldataout = output_data_gm;
                }
                // WM case
                if (o == 2)
                {
                    fsub = "_wm";
                    asldataout = output_data_wm;
                }
            }

            // Case to handle extrapolation output
            if (opts.extrapolate_option.value())
            {
                asldataout = output_data_extrapolated;
            }

            //output data. Use input volume as basis for output volume to ensure consistent metadata
            if (opts.out.set())
            {
                Matrix outmtx;
                volume4D<float> dataout = data * 0;
                ;
                stdform2data(asldataout, outmtx, outblocked, outpairs);
                dataout.setmatrix(outmtx, mask);
                save_volume4D(dataout, opts.out.value() + fsub);
            }

            //take mean at each TI
            if (opts.meanout.set())
            {
                timeanout(asldataout, mask, opts.meanout.value() + fsub, outpairs);
                /*cout << "Outputting ASL data mean at each TI" << endl;
	  Matrix meanti(ntis,nvox);
	  for (int n=0; n<ntis; n++) 
	  {
	  meanti.Row(n+1)=mean(asldata[n],1);
	  }
	  volume4D<float> meanout = data*0;;
	  meanout.setmatrix(meanti,mask);
	  save_volume4D(meanout,opts.meanout.value());*/
            }

            //split data into separate file for each TI
            if (opts.splitout.set())
            {
                splitout(asldataout, mask, opts.splitout.value() + fsub);
            }

            //do epochwise output
            if (opts.epochout.set())
            {
                int eplen = opts.epochlen.value();
                int epol = opts.epochover.value();

                bool tiunit; //indicates if we want epochs over TIs rather than over repeats
                string epunit = opts.epochunit.value();
                if (epunit.compare("rpt") == 0)
                    tiunit = false;
                else if (epunit.compare("tis") == 0)
                    tiunit = true;
                else
                    throw Exception("Unrecognised epoch unit");

                //if (ispairs && !opts.tcdiff.value() ) throw Exception("Epoch ceration is currently incompatible with data that contains pairs if differencing has not been performed");
                /* NOT CORRECT - the pairs are stores within the TI
      if (outpairs) {
	//if we are dealing with pairs then we have two measurements for every repeat
	// epoch parameters are (meant to be) in number of repeats
	// hence multiply by 2 in this scenario
	eplen*=2;
	epol*=2;
	}*/

                epochout(asldataout, mask, opts.epochout.value() + fsub, eplen, epol, outpairs, tiunit);
                /*cout << "Creating ASL data epochs" << endl;
      int eplen = opts.epochlen.value();
      int epol = opts.epochover.value();
      if (epol>=eplen) throw Exception("The epoch overlap may not exceed or equal the length of the epoch");
      int epadv = eplen-epol;

      //NOTE that ispairs indicates that in asldata the measurements come in pairs that we collected at the same time
      //int ncollect=1; //number of measurements to collect for a single repeat-TI
      //if (ispairs) ncollect=2;

      int e=1;
      Matrix epoch_temp(ntis,nvox);
      volume4D<float> epochout = data*0;;
      while (e*epadv<=nmeas)
	{
	  //Matrix ti_temp(eplen,nvox);
	  //go through the TIs
	  for (int ti=1; ti<=ntis; ti++) {
	    epoch_temp.Row(ti) = mean(asldata[ti-1].Rows((e-1)*epadv+1,e*epadv),1);
	  }

	  char cstr [5];
	  if (e<10) sprintf(cstr,"00%d",e);
	  else if (e<100) sprintf(cstr,"0%d",e);
	  else if (e<1000) sprintf(cstr,"%d",e);
	  else throw Exception("More than 1000 epochs in this ASL data file, sorry cannot handle this operation");
	  string epno(cstr);
	  
	  epochout.setmatrix(epoch_temp,mask);
	  save_volume4D(epochout,opts.epochout.value()+epno);

	  e++;
	}

      int unused=nmeas-(e-1)*epadv;
      cout << unused << endl;
      assert(unused>=0);
      if (unused>0) cout << "Number of measurements from end of data discarded: " << unused << endl;*/
            }
        }

        // OUTPUT: supplementary output section, things that are not compatible with the splitting of pairs into separate outputs

        // do deconvolution
        if (opts.deconvout.set())
        {
            //load aif
            volume4D<float> aif;
            read_volume4D(aif, opts.aif.value());
            Matrix aifmtx;
            aifmtx = aif.matrix(mask);

            deconvout(asldata, mask, aifmtx, opts.deconvout.value());
        }
    }

    catch (Exception &e)
    {
        cerr << endl
             << e.what() << endl;
        return 1;
    }
    catch (X_OptionError &e)
    {
        cerr << endl
             << e.what() << endl;
        return 1;
    }

    cout << "Done." << endl;

    return 0;
}
