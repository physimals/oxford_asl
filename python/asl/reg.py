#!/bin/env python

# ASL_REG: Registration for ASL data
#
# Michael Chappell, IBME QuBIc & FMRIB Image Analysis Groups
#
# Copyright (c) 2008-2016 University of Oxford
#
# SHCOPYRIGHT

import os, sys
import shutil
import tempfile
import collections
from optparse import OptionParser, OptionGroup

import asl.fslhelpers as fsl
import asl

def main():
    usage = """ASL_REG
    Registration for ASL data

    asl_reg -i <input> -s <struct> [options]"""

    p = OptionParser(usage=usage, version=asl.__version__)
    g = OptionGroup(p, "Required arguments")
    g.add_option("-i", dest="input", help="input image - e.g. perfusion-weighted image")
    g.add_option("-s", dest="struct", help="structural brain image - wholehead")
    p.add_option_group(g)

    g = OptionGroup(p, "Optional arguments")
    g.add_option("-o", dest="outdir", help="output directory", default=os.getcwd())
    g.add_option("--sbet", dest="sbet", help="structural brain image - brain extracted")
    g.add_option("--init", dest="init", help="initial transformation matrix for input to structural image")
    p.add_option_group(g)

    g = OptionGroup(p, "Extra 'final' registration refinement")
    g.add_option("-c", dest="cfile", help="ASL control/calibration image for initial registration - brain extracted")
    g.add_option("-m", dest="mask", help="brain mask for brain extraction of the input image")
    g.add_option("--tissseg", dest="tissseg", help="tissue segmenation image for bbr (in structural image space)")
    g.add_option("--finalonly", action="store_true", dest="finalonly", help="only run the 'final' registration step", default=False)
    p.add_option_group(g)

    g = OptionGroup(p, "Distortion correction using fieldmap (see epi_reg)")
    g.add_option("--fmap", dest="fmap", help="fieldmap image (in rad/s)")
    g.add_option("--fmapmag", dest="fmapmag", help="fieldmap magnitude image - wholehead extracted")
    g.add_option("--fmapmagbrain", dest="fmapmagbrain", help="fieldmap magnitude image - brain extracted")
    g.add_option("--wmseg", dest="wmseg", help="white matter segmentation of T1 image")
    g.add_option("--echospacing", dest="echospacing", help="Effective EPI echo spacing (sometimes called dwell time) - in seconds", type="float")
    g.add_option("--pedir", dest="pedir", help="phase encoding direction, dir = x/y/z/-x/-y/-z")
    g.add_option("--nofmapreg", dest="nofmapreg", help="do not perform registration of fmap to T1 (use if fmap already registered)", action="store_true", default=False)
    p.add_option_group(g)

    g = OptionGroup(p, "Deprecated")
    g.add_option("-r", dest="lowstruc", help="extra low resolution structural image - brain extracted")
    g.add_option("--inweight", dest="inweight", help="specify weights for input image - same functionality as the flirt -inweight option", type="float")
    p.add_option_group(g)

    options, args = p.parse_args()

    print("ASL_REG")

    vars = collections.defaultdict(str)
    vars["fsldir"] = os.environ["FSLDIR"]

    vars["infile"] = fsl.Image(options.input, "input").full
    print("Input file is: %(infile)s" % vars)

    vars["outdir"] = options.outdir
    print("Output directory is: %(outdir)s" % vars)
    
    if not os.path.exists(vars["outdir"]):
        print("Creating output directory")
        fsl.mkdir(vars["outdir"])
  
    if options.debug:
        vars["tempdir"] = os.path.join(os.getcwd(), "tmp_asl_reg")
        fsl.mkdir(vars["tempdir"])
    else:
        vars["tempdir"] = tempfile.mkdtemp("_asl_reg")
    print("Using temporary dir: %(tempdir)s" % vars)

    # Initial matrix option
    if options.init:
        vars["init"] = "-init %s" % options.init
        vars["epi_init"] = "--init=%s" % options.init

    # Weighting applied to input image
    if options.inweight:
        vars["inweight"] = "-inweight %s" % options.inweight

    # Degrees of freedom - we will routinely use 6
    vars["dof"] = 6
    vars["xyssch"] = os.path.join(vars["fsldir"], "etc", "flirtsch", "xyztrans.sch")

    # Optional flirt schedule for main transformation of asl to structural
    if options.flirtsch:
        print("Using supplied FLIRT schedule")
        vars["flirtsch"] = options.flirtsch
    else:
        vars["flirtsch"] = os.path.join(vars["fsldir"], "etc", "flirtsch", "simple3D.sch")

    # BET the structural image if required
    if options.sbet:
        vars["sbet"] = fsl.Image("options.sbet", "sbet").full
    else:
        print("Running BET on structural image")
        fsl.bet("%(struc)s %(tempdir)s/struc_brain" % vars)
        vars["sbet"] = "%(tempdir)s/struc_brain" % vars

    # do the MAIN registration run - use the supplementary image for this if available
    if not options.finalonly:
        print("Registration MAIN stage (FLIRT)")

        # Check if a supplementary image has been provded on which to base (inital) registration
        if options.cfile:
            vars["cfile"] = options.cfile
        else:
            vars["cfile"] = vars["infile"]
            
        if not options.lowstruc:
            # Step1: 3DOF translation only transformation
            fsl.flirt("-in %(cfile)s -ref %(sbet)s -schedule %(xyzsch)s -omat %(tempdir)s/low2high1.mat -out %(tempdir)s/low2hig1 %(init)s %(inweight)s" % vars)
            # Step2: 6DOF transformation with small search region
            fsl.flirt("-in %(cfile)s -ref %(sbet)s -dof %(dof)s -omat %(tempdir)s/low2high.mat -init %(tempdir)s/low2high1.mat -schedule %(flirtsch)s -out %(tempdir)s/low2high  %(inweight)s" % vars)
        else:
            # We have a structural image in perfusion space use it to improve registration
            vars["lowstruc"] = fsl.Image(options.lowstruc, "lowstruc").full
            print("Using structral image in perfusion space (%(lowstruc)s)" % vars)
            # Step1: 3DOF translation only transformation perfusion->lowstruc
            fsl.flirt("-in %(cfile)s -ref %(lowstruc)s -schedule %(xyzsch)s -omat %(tempdir)s/low2low1.mat %(init)s %(inweight)s" % vars)
            # Step2: 6DOF limited transformation in perfusion space
            fsl.flirt("-in %(cfile)s -ref %(lowstruc)s -dof %(dof)s -schedule %(flirtsch)s -init %(tempdir)s/low2low1.mat -omat %(tempdir)s/low2low.mat %(inweight)s" % vars)
            # Step3: 6DOF transformation of lowstruc to struc
            fsl.flirt("-in %(lowstruc)s -ref %(sbet)s -omat %(tempdir)s/str2str.mat" % vars)
            # Step4: combine the two transformations
            fsl.Prog("convert_xfm")("-omat %(tempdir)s/low2high.mat -concat %(tempdir)s/str2str.mat %(tempdir)s/low2low.mat" % vars)

        # update the init text ready for the 'final' step to start with the result of the MAIN registration
        vars["epi_init"] = "--init=%(tempdir)s/low2high.mat" % vars

        # OUTPUT from MAIN registration
        shutl.copy("%(tempdir)s/low2high.mat" % vars, "%(outdir)s/asl2struct.mat" % vars)

    # do the FINAL registration run using BBR - this ONLY makes sense when the input is a perfusion image (or something with decent tissue contrast)
    if not options.mainonly:
        print("Registration FINAL stage (BBR)")

        if options.tissseg:
            vars["tissseg"] = options.tissseg
        else:
            # Running FAST segmentation
            fsl.fast("-o %(tempdir)s/struct_fast %(sbet)s" % vars)
            # WM segmentation
            fsl.maths("%(tempdir)s/struct_fast_pve_2 -thr 0.5 -bin ${tempdir}/fast_tissseg" % vars)
            vars["tissseg"] = "${tempdir}/fast_tissseg" % vars

        # Brain extract the perfusion image - using supplied mask or mask derived from the strctural BET
        if options.mask:
            vars["mask"] = options.mask
        else:
            fsl.Prog("convert_xfm")("-omat %(tempdir)s/high2low.mat -inverse %(tempdir)s/low2high.mat" % vars)

            fsl.maths("%(sbet)s -thr 0 -bin %(tempdir)s/struct_brain_mask" % vars)
            fsl.flirt("-in %(tempdir)s/struct_brain_mask -ref %(infile)s -applyxfm -init %(tempdir)s/high2low.mat -out %(tempdir)s/mask -interp trilinear" % vars)
            fsl.maths("%(tempdir)s/mask -thr 0.25 -bin -fillh %(tempdir)s/mask" % vars)
            fsl.Prog("fslcpgeom")("%(infile)s %(tempdir)s/mask" % vars)
            vars["mask"] = "%(tempdir)s/mask" % vars

        # Apply mask to asldata
        fsl.maths("%(infile)s -mas %(mask)s %(tempdir)s/asldata_brain" % vars)

        # Copy mask to output for future reference
        fsl.imcp("%(mask)s" % vars, "%(outdir)s/mask" % vars)
        
        # Do a final refinement of the registration using the perfusion and the white matter 
        # segmentation - using epi_reg to get BBR (and allow for fieldmap correction in future)
        if options.fmap:
            if options.nofmapreg:
                vars["fmapregstr"] = "--nofmapreg"
            fsl.Prog("epi_reg")("--epi=%(tempdir)s/asldata_brain --t1=%(struc)s --t1brain=%(sbet)s $epi_inittext --out=%(tempdir)s/low2high_final --wmseg=%(tissseg)s %(inweight)s --fmap=%(fmap)s --fmapmag=%(fmapmag)s --fmapmagbrain=%(fmapmagbrain)s --pedir=%(pedir)s --echospacing=%(echospacing)s %(fmapregstr)s" % vars)
        else:
            fsl.Prog("epi_reg")("--epi=%(tempdir)s/asldata_brain --t1=%(struc)s --t1brain=%(sbet)s $epi_inittext --out=%(tempdir)s/low2high_final --wmseg=%(tissseg)s %(inweight)s" % vars)
                
        #flirt -ref %(sbet)s -in %(infile)s -dof 6 -cost bbr -wmseg $wmseg -init %(tempdir)s/low2high.mat -omat %(tempdir)s/low2high.mat -out %(tempdir)s/low2high_final -schedule ${FSLDIR}/etc/flirtsch/bbr.sch
        print("BBR end")

        print("Saving FINAL output")
        if not options.finalonly:
            # Save the initial transformation matrix to allow chekcing if this part failed
            shutil.copy("%(outdir)s/asl2struct.mat" % vars, "%(outdir)s/asl2struct_init.mat" % vars)

        # The transformation matrix from epi_reg - this overwrites the version from MAIN registration
        shutil.copy("%(tempdir)s/low2high_final.mat" % vars, "%(outdir)s/asl2struct.mat" % vars)

        # Often useful to have the inverse transform, so calcuate it
        fsl.Prog("convert_xfm")("-omat %(outdir)s/struct2asl.mat -inverse %(outdir)s/asl2struct.mat" % vars)

        if options.fmap:
            # The warp from epi_reg
            fsl.imcp("%(tempdir)s/low2high_final_warp" % vars, "%(outdir)s/asl2struct_warp" % vars)

        # Save the transformed image to check on the registration
        fsl.imcp("%(tempdir)s/low2high_final" % vars, "%(outdir)s/asl2struct" % vars)
        
        # Copy the edge image from epi_reg output as that is good for visualisation
        fsl.imcp("%(tissseg)s" % vars, "%(outdir)s/tissseg" % vars)
        fsl.imcp("%(tempdir)s/low2high_final_fast_wmedge" % vars, "%(outdir)s/tissedge" % vars)

    # ASL-->standard transformation (if specified)
    #if options.transflag
    #    print("Combining transformations")
    #    fsl.Prog("convert_xfm")("-omat %(outdir)s/asl2std.mat -concat $trans %(tempdir)s/low2high.mat" % vars)

    #if options.lowstruc:
    #    ASL--> low structral transformtaion (if supllied)
    #    shutil.copy("%(tempdir)s/low2low.mat" % vars, "%(outdir)s/asl2lowstruct.mat" % vars)

    print("ASL_REG - Done.")

if __name__ == "__main__":
    main()

