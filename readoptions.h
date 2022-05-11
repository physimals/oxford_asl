/*  readoptions.h

    Michael Chappell - FMRIB Image Analysis Group

    Moss Zhao - IBME Quantitative Biomedical Inference (QuBIc) Group

    Copyright (C) 2015 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(ReadOptions_h)
#define ReadOptions_h

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "utils/log.h"
#include "utils/options.h"

namespace OXASL
{
class ReadOptions
{
public:
    static ReadOptions &getInstance();
    ~ReadOptions() { delete ropt; }
    Utilities::Option<bool> help;
    Utilities::Option<bool> version;

    Utilities::Option<std::string> datafile;
    Utilities::Option<std::string> maskfile;
    Utilities::Option<int>         ntis;

    Utilities::Option<std::string> inblockform;
    Utilities::Option<std::string> inaslform;
    Utilities::Option<std::string> rpts;
    Utilities::Option<bool>        ispairs;

    Utilities::Option<bool> splitpairs;
    Utilities::Option<bool> tcdiff;
    Utilities::Option<bool> surrtcdiff;

    Utilities::Option<std::string> outblockform;
    Utilities::Option<std::string> out;

    Utilities::Option<std::string> meanout;
    Utilities::Option<std::string> splitout;

    Utilities::Option<std::string> epochout;
    Utilities::Option<int>         epochlen;
    Utilities::Option<int>         epochover;
    Utilities::Option<std::string> epochunit;

    Utilities::Option<std::string> deconvout;
    Utilities::Option<std::string> aif;

    // Extrapolate the edge of the brain to fix the artefact on the edge of the brain
    // Assumes an eroded brain
    Utilities::Option<bool> extrapolate_option;
    Utilities::Option<int>  neighbour;

    // Partial volume correction (linear regression method) parameters
    Utilities::Option<std::string> pv_gm_file;
    Utilities::Option<std::string> pv_wm_file;
    Utilities::Option<int>         kernel;

    bool parse_command_line(int argc, char **argv);

private:
    ReadOptions();
    const ReadOptions &operator=(ReadOptions &);
    ReadOptions(ReadOptions &);

    Utilities::OptionParser options;

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

    help(std::string("-h,--help"), false,
        std::string("display this message"),
        false, Utilities::no_argument)
    ,

    version(std::string("-v,--version"), false,
        std::string("display version identification"),
        false, Utilities::no_argument)
    ,

    //input files
    datafile(std::string("--data,--datafile"), std::string("ASL datafile"),
        std::string("data file"),
        true, Utilities::requires_argument)
    , maskfile(std::string("--mask"), std::string("maskfile\n"),
          std::string("mask"),
          false, Utilities::requires_argument)
    ,

    //input file information
    ntis(std::string("--ntis"), 0,
        std::string("Number of TIs in file"),
        true, Utilities::requires_argument)
    , inblockform(std::string("--ibf,--inblockform"), std::string("rpt"),
          std::string("Input block format:\n          rpt - blocks of measurements that include all TIs\n          tis - blocks of repeated measurements at a single TI"),
          false, Utilities::requires_argument)
    , inaslform(std::string("--iaf,--inaslform"), std::string("diff"),
          std::string("ASL data form:\n          diff - differenced data {default}\n          tc   - Tag-Control pairs\n          ct   - Control-Tag pairs\n          tcb  - Tag-Control pairs, tags and controls grouped together within block\n          ctb - Control-Tag pairs, tags and controls grouped together within block"),
          false, Utilities::requires_argument)
    , rpts(std::string("--rpts"), std::string("NULL"),
          std::string("Number of repeats at each TI as comma separated list, not required if the number of repeats is same for all TIs  (only for use with --ibf=tis)"),
          false, Utilities::requires_argument)
    ,

    ispairs(std::string("--pairs,--inpairs"), false,
        std::string("Data contains adjacent pairs of measuremnts (e.g. Tag, Control) DEPRECEATED used --iaf instead\n"),
        false, Utilities::no_argument)
    ,

    //asaq(std::string("--asaq"),false,
    //	std::string("Data is as aquired: same as --blocked --pairs"),
    //	false,Utilities::no_argument),

    // manipulation options
    splitpairs(std::string("--spairs"), false,
        std::string("Split the pairs within the data, e.g. to separate tag and control images in output"),
        false, Utilities::no_argument)
    , tcdiff(std::string("--diff"), false,
          std::string("Take the difference between the pairs, i.e. Tag control difference"),
          false, Utilities::no_argument)
    , surrtcdiff(std::string("--surrdiff"), false,
          std::string("Do surround subtraction on the pairs\n"),
          false, Utilities::no_argument)
    ,

    //basic output
    outblockform(std::string("--obf,--outblockform"), std::string("notset"),
        std::string("Output block format (for --out=):\n          rpt - blocks of measurements that include all TIs\n          tis - blocks of repeated measurements at a single TI\n          Default is same as input block format (--ibf)"),
        false, Utilities::requires_argument)
    , out(std::string("--out"), std::string("Out filename"),
          std::string("Output data file"),
          false, Utilities::requires_argument)
    ,

    // other output options
    meanout(std::string("--mean"), std::string(""),
        std::string("Output ASL data having taken mean at each TI to file"),
        false, Utilities::requires_argument)
    , splitout(std::string("--split"), std::string(""),
          std::string("Split data into separate files each each TI, specify filename root\n"),
          false, Utilities::requires_argument)
    ,

    epochout(std::string("--epoch"), std::string(""),
        std::string("Output epochs of ASL data (takes mean at each TI within the epoch)"),
        false, Utilities::requires_argument)
    , epochlen(std::string("--elen,--epochlen"), 1,
          std::string("Length of epochs in number of repeats"),
          false, Utilities::requires_argument)
    , epochover(std::string("--eol,--epochol"), 0,
          std::string("Ammount of overlap between epochs in number of repeats"),
          false, Utilities::requires_argument)
    , epochunit(std::string("--eunit,--epochunit"), std::string("rpt"),
          std::string("Epochs to be determined over:\n          rpt - repeats in the data {default}\n          tis - TIs in the data\n"),
          false, Utilities::requires_argument)
    ,

    deconvout(std::string("--deconv"), std::string(""),
        std::string("Deconvolution of data with arterial input functions"),
        false, Utilities::requires_argument)
    , aif(std::string("--aif"), std::string(""),
          std::string("Arterial input functions for deconvolution (4D volume, one aif for each voxel within mask)\n"),
          false, Utilities::requires_argument)
    ,

    // Extrapolate the edge of the brain to fix the artefact on the edge of the brain
    // Assumes an eroded brain
    extrapolate_option(std::string("--extrapolate"), false,
        std::string("Option to extrapolate the edge of the brain to fix the artefact on the edge of the brain"),
        false, Utilities::no_argument)
    , neighbour(std::string("--neighbour"), 5, std::string("Neighbour size for extrapolation, must be an odd number between 3 and 9. Default: 5\n"),
          false, Utilities::requires_argument)
    ,

    // Partial volume (linear regression) options
    pv_gm_file(std::string("--pvgm"), std::string(""), std::string("GM partial volume map"),
        false, Utilities::requires_argument)
    , pv_wm_file(std::string("--pvwm"), std::string(""), std::string("WM partial volume map"),
          false, Utilities::requires_argument)
    , kernel(std::string("--kernel"), 5, std::string("Kernel size (in voxels) of partial volume correction, must be an odd number between 3 and 9. Default: 5\n"),
          false, Utilities::requires_argument)
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

    catch (Utilities::X_OptionError &e)
    {
        options.usage();
        std::cerr << std::endl
             << e.what() << std::endl;
    }
    catch (std::exception &e)
    {
        std::cerr << e.what() << std::endl;
    }
}
}
#endif
