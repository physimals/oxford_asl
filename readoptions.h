/*  readoptions.h

    Michael Chappell - FMRIB Image Analysis Group

    Moss Zhao - IBME Quantitative Biomedical Inference (QuBIc) Group

    Copyright (C) 2015 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(ReadOptions_h)
#define ReadOptions_h

#include "utils/log.h"
#include "utils/options.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace Utilities;

namespace OXASL
{
class ReadOptions
{
public:
    static ReadOptions &getInstance();
    ~ReadOptions() { delete ropt; }
    Option<bool> help;
    Option<bool> version;

    Option<string> datafile;
    Option<string> maskfile;
    Option<int> ntis;

    Option<string> inblockform;
    Option<string> inaslform;
    Option<string> rpts;
    Option<bool> ispairs;

    Option<bool> splitpairs;
    Option<bool> tcdiff;
    Option<bool> surrtcdiff;

    Option<string> outblockform;
    Option<string> out;

    Option<string> meanout;
    Option<string> splitout;

    Option<string> epochout;
    Option<int> epochlen;
    Option<int> epochover;
    Option<string> epochunit;

    Option<string> deconvout;
    Option<string> aif;

    // Extrapolate the edge of the brain to fix the artefact on the edge of the brain
    // Assumes an eroded brain
    Option<bool> extrapolate_option;
    Option<int> neighbour;

    // Partial volume correction (linear regression method) parameters
    Option<string> pv_gm_file;
    Option<string> pv_wm_file;
    Option<int> kernel;

    bool parse_command_line(int argc, char **argv);

private:
    ReadOptions();
    const ReadOptions &operator=(ReadOptions &);
    ReadOptions(ReadOptions &);

    OptionParser options;

    static ReadOptions *ropt;
};

inline ReadOptions &ReadOptions::getInstance()
{
    if (ropt == NULL)
        ropt = new ReadOptions();

    return *ropt;
}

inline ReadOptions::ReadOptions()
    :

    help(string("-h,--help"), false,
        string("display this message"),
        false, no_argument)
    ,

    version(string("-v,--version"), false,
        string("display version identification"),
        false, no_argument)
    ,

    //input files
    datafile(string("--data,--datafile"), string("ASL datafile"),
        string("data file"),
        true, requires_argument)
    , maskfile(string("--mask"), string("maskfile\n"),
          string("mask"),
          false, requires_argument)
    ,

    //input file information
    ntis(string("--ntis"), 0,
        string("Number of TIs in file"),
        true, requires_argument)
    , inblockform(string("--ibf,--inblockform"), string("rpt"),
          string("Input block format:\n          rpt - blocks of measurements that include all TIs\n          tis - blocks of repeated measurements at a single TI"),
          false, requires_argument)
    , inaslform(string("--iaf,--inaslform"), string("diff"),
          string("ASL data form:\n          diff - differenced data {default}\n          tc   - Tag-Control pairs\n          ct   - Control-Tag pairs\n          tcb  - Tag-Control pairs, tags and controls grouped together within block\n          ctb - Control-Tag pairs, tags and controls grouped together within block"),
          false, requires_argument)
    , rpts(string("--rpts"), string("NULL"),
          string("Number of repeats at each TI as comma separated list, not required if the number of repeats is same for all TIs  (only for use with --ibf=tis)"),
          false, requires_argument)
    ,

    ispairs(string("--pairs,--inpairs"), false,
        string("Data contains adjacent pairs of measuremnts (e.g. Tag, Control) DEPRECEATED used --iaf instead\n"),
        false, no_argument)
    ,

    //asaq(string("--asaq"),false,
    //	string("Data is as aquired: same as --blocked --pairs"),
    //	false,no_argument),

    // manipulation options
    splitpairs(string("--spairs"), false,
        string("Split the pairs within the data, e.g. to separate tag and control images in output"),
        false, no_argument)
    , tcdiff(string("--diff"), false,
          string("Take the difference between the pairs, i.e. Tag control difference"),
          false, no_argument)
    , surrtcdiff(string("--surrdiff"), false,
          string("Do surround subtraction on the pairs\n"),
          false, no_argument)
    ,

    //basic output
    outblockform(string("--obf,--outblockform"), string("notset"),
        string("Output block format (for --out=):\n          rpt - blocks of measurements that include all TIs\n          tis - blocks of repeated measurements at a single TI\n          Default is same as input block format (--ibf)"),
        false, requires_argument)
    , out(string("--out"), string("Out filename"),
          string("Output data file"),
          false, requires_argument)
    ,

    // other output options
    meanout(string("--mean"), string(""),
        string("Output ASL data having taken mean at each TI to file"),
        false, requires_argument)
    , splitout(string("--split"), string(""),
          string("Split data into separate files each each TI, specify filename root\n"),
          false, requires_argument)
    ,

    epochout(string("--epoch"), string(""),
        string("Output epochs of ASL data (takes mean at each TI within the epoch)"),
        false, requires_argument)
    , epochlen(string("--elen,--epochlen"), 1,
          string("Length of epochs in number of repeats"),
          false, requires_argument)
    , epochover(string("--eol,--epochol"), 0,
          string("Ammount of overlap between epochs in number of repeats"),
          false, requires_argument)
    , epochunit(string("--eunit,--epochunit"), string("rpt"),
          string("Epochs to be determined over:\n          rpt - repeats in the data {default}\n          tis - TIs in the data\n"),
          false, requires_argument)
    ,

    deconvout(string("--deconv"), string(""),
        string("Deconvolution of data with arterial input functions"),
        false, requires_argument)
    , aif(string("--aif"), string(""),
          string("Arterial input functions for deconvolution (4D volume, one aif for each voxel within mask)\n"),
          false, requires_argument)
    ,

    // Extrapolate the edge of the brain to fix the artefact on the edge of the brain
    // Assumes an eroded brain
    extrapolate_option(string("--extrapolate"), false,
        string("Option to extrapolate the edge of the brain to fix the artefact on the edge of the brain"),
        false, no_argument)
    , neighbour(string("--neighbour"), 5, string("Neighbour size for extrapolation, must be an odd number between 3 and 9. Default: 5\n"),
          false, requires_argument)
    ,

    // Partial volume (linear regression) options
    pv_gm_file(string("--pvgm"), string(""), string("GM partial volume map"),
        false, requires_argument)
    , pv_wm_file(string("--pvwm"), string(""), string("WM partial volume map"),
          false, requires_argument)
    , kernel(string("--kernel"), 5, string("Kernel size (in voxels) of partial volume correction, must be an odd number between 3 and 9. Default: 5\n"),
          false, requires_argument)
    ,

    options("asl_file", "asl_file --data=<asldata> --ibf=rpt --iaf=tc --diff --out=<diffdata>\n")
{
    try
    {
        options.add(help);
        options.add(version);

        options.add(datafile);
        options.add(maskfile);

        options.add(ntis);
        options.add(inblockform);
        options.add(inaslform);
        options.add(rpts);
        options.add(ispairs);

        options.add(splitpairs);
        options.add(tcdiff);
        options.add(surrtcdiff);

        options.add(extrapolate_option);
        options.add(neighbour);

        options.add(pv_gm_file);
        options.add(pv_wm_file);
        options.add(kernel);

        options.add(out);
        options.add(outblockform);

        options.add(meanout);
        options.add(splitout);

        options.add(epochout);
        options.add(epochlen);
        options.add(epochover);
        options.add(epochunit);

        options.add(deconvout);
        options.add(aif);
    }

    catch (X_OptionError &e)
    {
        options.usage();
        cerr << endl
             << e.what() << endl;
    }
    catch (std::exception &e)
    {
        cerr << e.what() << endl;
    }
}
}
#endif
