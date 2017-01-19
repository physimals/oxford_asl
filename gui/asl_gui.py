#!/usr/bin/env python
#

from collections import OrderedDict

import os

import wx
import props

def conditional(p1, p2):
    return props.Widget(p1, visibleWhen=lambda i: getattr(i, p2))

class AslOptions(props.HasProperties):

    # Data tab
    inputImage = props.FilePath(exists=True, required=True)
    inversionTimes = props.String()
    bolusDuration = props.Real(default=1, minval=0, maxval=100)
    labelling = props.Choice(["pASL", "cASL/pcASL"])
    tagControlPairs = props.Boolean(default=True)
    controlFirst = props.Boolean(default=False)
    dataOrder = props.Choice(["Repeats", "TIs"])
    staticTissue = props.Choice(["normal", "background suppressed", "pre-saturation"]) # Only if tagControlPairs
    useStructuralImage = props.Boolean(default=False)
    structuralImage = props.FilePath(exists=True, required=False) # Only if useStructuralImage
    runBet = props.Boolean(default=True) # Only if useStructuralImage

    # Analysis tab
    outputDir = props.FilePath(exists=False, required=True, isFile=False)
    brainMask = props.FilePath(exists=True, required=False)
    paramVariance = props.Boolean(default=False)
    bolusArrivalTime = props.Real(default=0.7, minval=0, maxval=100)
    t1 = props.Real(default=1.3, minval=0, maxval=100)
    t1b = props.Real(default=1.6, minval=0, maxval=100)
    inversionEffic = props.Real(default=0.98, minval=0, maxval=100)
    spatialSmoothing = props.Boolean(default=False)
    inferT1 = props.Boolean(default=False)
    includeVascular = props.Boolean(default=False)
    fixBolusDuration = props.Boolean(default=True)

    # Registration tab
    useSpaceTransform = props.Boolean(default=False) # Only if useStructuralImage
    spaceTransform = props.FilePath(exists=True, required=False)  # Only if useSpaceTransform
    useAlternateStandardImage = props.Boolean(default=False)  # Only if useStructuralImage
    alternateStandardImage = props.FilePath(exists=True, required=False)  # Only if useAlternateStandardImage

    # Calibration tab
    doCalibration = props.Boolean(default=False)
    calibrateMode = props.Choice(["Long TR", "Saturation Recovery"])
    m0Image = props.FilePath(exists=True, required=False)  # Only if useCoilSensitivityRefImage
    useCoilSensitivityRefImage = props.Boolean(default=False)
    coilSensitivityRefImage = props.FilePath(exists=True, required=False) # Only if useCoilSensitivityRefImage
    calibrateGain = props.Real(default=1.0, minval=0, maxval=100)
    refTissueType = props.Choice(["CSF", "White matter", "Grey matter", "none"])
    useRefTissueMask = props.Boolean(default=True) #?? Only if refTissueType != none???
    refTissueMask = props.FilePath(exists=True, required=False)
    refT1 = props.Real(default=4.3, minval=0, maxval=100)
    refT2 = props.Real(default=0.75, minval=0, maxval=100)
    bloodT2 = props.Real(default=0.15, minval=0, maxval=100)
    seqTR = props.Real(default=3.2, minval=0, maxval=100)
    seqTE = props.Real(default=0.0, minval=0, maxval=100)

    #def setOutputImage(self, value, valid, *a):
    #    if not valid: return
    #    value = removeExt(value, ['.nii.gz', '.nii'])
    #    self.outputImage = value + '_brain'

    def writableDir(self, d):
        if os.path.isdir(d):
            if not os.access(d, os.R_OK | os.W_OK | os.X_OK):
                raise RuntimeError('Output dir exists but is not accessible')
        else:
            os.mkdir(d)

    def fslProg(self, prog):
        return os.path.join(os.environ["FSLDIR"], "bin/%s" % prog)

    def imtest(self, f):
        """
        Test if file contains an image using external FSL tool. Not
        great and native solution would be better
        """
        ret = os.system("%s %s" % (self.fslProg("imtest"), f))
        return ret == 1

    def runcmd(self, c):
        """
        Run a command
        """
        print(c)
        ret = os.system(c)
        if ret != 0:
            raise RuntimeError("Error executing command:\n%s" % c)

    def warn(self, message, caption="Warning"):
        dlg = wx.MessageDialog(None, message, caption, wx.OK | wx.ICON_WARNING)
        dlg.ShowModal()
        dlg.Destroy()

    def runAsl(self):
        try:
            errors = self.validateAll()
            if len(errors) > 0:
                errs = ["{}: {}".format(prop, msg) for prop, msg in errors]
                raise RuntimeError("\n".join(errs))

            # Make output dir FIXME what if already exists? Currently overwrite
            self.writableDir(self.outputDir)
            nativeSpaceDir = os.path.join(self.outputDir, "native_space")
            self.writableDir(nativeSpaceDir)

            # Inversion times is a list, FIXME comma or space separated? FIXME check all numbers
            if self.inversionTimes is None or self.inversionTimes.strip() == "":
                raise RuntimeError("No inversion times specified")

            # FIXME check num tis against input image?
            tis = self.inversionTimes.replace(",", " ").split()
            print(tis)

            aslFileCommand = [self.fslProg("asl_file"),]
            aslFileCommand.append("--data=%s" % self.inputImage)
            aslFileCommand.append("--out=%s/diffData" % nativeSpaceDir)
            aslFileCommand.append("--obf=rpt --ntis=%i" % len(tis))

            if self.dataOrder == "Repeats":
                aslFileCommand.append("--ibf=rpt")
            else:
                aslFileCommand.append("--ibf=tis")

            if self.controlFirst:
                aslDataForm = "--iaf=ct"
            else:
                aslDataForm = "--iaf=tc"

            # 0 "none" 1 "pairwise"
            if self.tagControlPairs:
                aslFileCommand.append("%s --diff" % aslDataForm)
            else:
                aslFileCommand.append("--iaf=diff")

            self.runcmd(" ".join(aslFileCommand))

            aslCommand = [self.fslProg("oxford_asl"),]
            aslCommand.append("-i %s/diffData" % nativeSpaceDir)
            aslCommand.append("-o %s" % self.outputDir)
            aslCommand.append("--tis %s" % self.inversionTimes)
            if self.paramVariance: aslCommand.append("--vars")
            aslCommand.append("--bolus %f" % self.bolusDuration)
            aslCommand.append("--bat %f" % self.bolusArrivalTime)
            aslCommand.append("--t1 %f" % self.t1)
            aslCommand.append("--t1b %f" % self.t1b)
            aslCommand.append("--alpha %f" % self.inversionEffic)
            if self.spatialSmoothing: aslCommand.append("--spatial")
            if self.inferT1: aslCommand.append("--infert1")
            if not self.includeVascular: aslCommand.append("--artoff")
            if self.fixBolusDuration: aslCommand.append("--fixbolus")

            if self.brainMask: aslCommand.append("-m %s" % self.brainMask)
            if self.labelling != "pASL": aslCommand.append("--casl")

            if self.useStructuralImage:
                if runBet:
                    strucProg = self.fslProg("bet")
                else:
                    strucProg = self.fslProg("imcp")
                strucCmd = "%s %s %s/structural_brain" % (strucProg, self.structuralImage, self.outputDir)
                self.runcmd(strucCmd)
                aslCommand.append("-s %d/structural_brain" % self.outputDir)

                if self.useAlternateStandardImage:
                    aslCommand.append("-t %s -S %s" % (self.spaceTransform, self.alternateStandardImage))

            foundRegTarget = False
            if self.tagControlPairs and self.staticTissue in  ("pre-saturation", "background suppressed"):
                aslFileCmd = [self.fslProg("asl_file"),]
                aslFileCmd.append("--data=%s" % self.inputImage)
                aslFileCmd.append("--ntis=%i" % len(tis))
                aslFileCmd.append("--spairs")
                aslFileCmd.append(aslDataForm)
                aslFileCmd.append("--out=%s/asldata_mc" % self.outputDir)
                self.runcmd(" ".join(aslFileCmd))

                mode = self.tagControlPairs and self.staticTissue == "normal"

                if self.staticTissue == "pre-saturation":
                    self.runcmd("%s %s/asldata_mc_odd %s/asldata_mc_odd_brain" % (self.fslProg("bet"), self.outputDir, self.outputDir))
                    aslCommand.append("--regfrom %s/asldata_mc_odd_brain" % self.outputDir)
                else:
                    self.runcmd("%s %s/asldata_mc_even %s/asldata_mc_even_brain" % (self.fslProg("bet"), self.outputDir, self.outputDir))
                    aslCommand.append("--regfrom %s/asldata_mc_even_brain" % self.outputDir)
                foundRegTarget = True


            if self.doCalibration and not foundRegTarget:
                self.runcmd("%s %s %s/MOimage_brain" % (
                self.fslProg("bet"), self.m0Image, self.outputDir))
                aslCommand.append("--regfrom  %s/MOimage_brain" % self.outputDir)

            self.runcmd(" ".join(aslCommand))

            if self.doCalibration:
                aslCalibCommand = [self.fslProg("asl_calib"),]
                aslCalibCommmand.append("-i %s/native_space/perfusion" % self.outputDir)
                aslCalibCommand.append("--tissref %s" % self.tissueType)
                aslCalibCommand.append("--t1r %f" % self.referenceT1)
                aslCalibCommand.append("--t2r %s " % self.referenceT2)
                aslCalibCommand.append("--t2b %f" % self.bloodT2)
                aslCalibCommand.append("--te %f" % self.seqTE)
                aslCalibCommand.append("-o %s/calibration" % self.outputDir)

                if self.useRefTissueMask: aslCalibCommand.append("-m %s" % self.refTissueMask)
                else:
                    aslCalibCommand.append("-s %s/structural_brain" % self.outputDir)
                    aslCalibCommand.append("-t %s/native_space/asl2struct.mat" % self.outputDir)

                if self.brainMask:
                    aslCalibCommand.append("--bmask %s" % self.brainMask)
                if self.calibrateMode == 0: # long TR
                    calibImage = self.m0Image
                    if self.tagControlPairs and staticTissue == "normal":
                        calibImage = self.outputDir + "/asldata_mc_even"
                    aslCalibCommand.append("-c %s" % calibImage)
                    aslCalibCommand.append("--mode longtr --tr %f" % self.seqTR)
                    aslCalibCommand.append("--cgain %f" % self.calibrateGain)
                    if self.useCoil: aslCalibCommand.append("--cref %s" % self.coilSensitivityRefImage)
                elif self.calibrateMode == 1: # Saturation Recovery - maybe change aslcalib -c option to even - for MC to think about!!!
                    aslCalibCommand.append("--tis %s" % self.inversionTimes)
                    aslCalibCommand.append("-c %s/asldata_mc_odd" % self.outputDir)

                    self.runcmd(" ".join(aslCalibCommand))

            if self.useCoilSensitivityRefImage:
                aslDivCommand = "-div %s" % self.coilSensitivityRefImage
            else:
                aslDivCommand = ""

            if self.imtest(self.outputDir + "/perfusion"):
                M0 = self.fileContents("%s/calibration/M0.txt" % self.outputDir)
                self.runcmd("%s %s/perfusion -div %s %s -mul 6000 %s/perfusion_calib" % (self.fslProg("fslmaths"), M0, self.aslDivCommand, self.outputDir))

            if self.imtest(self.outputDir + "/standard_space/perfusion"):
                M0 = self.fileContents("%s/calibration/M0.txt" % self.outputDir)
                self.runcmd("%s %s/standard_space/perfusion -div %s %s -mul 6000 %s/standard_space/perfusion_calib" % (self.fslProg("fslmaths"), M0, self.aslDivCommand, self.outputDir))

            if self.imtest(self.outputDir + "/structural_space/perfusion"):
                M0 = self.fileContents("%s/calibration/M0.txt" % self.outputDir)
                self.runcmd("%s %s/structural_space/perfusion -div %s %s -mul 6000 %s/structural_space/perfusion_calib" % (self.fslProg("fslmaths"), M0, self.aslDivCommand, self.outputDir))

            if self.useAlternateStandardImage:
                regOptions = "--prefix=%s/native_space/asl2struct.mat -r %s " % (self.outputDir, self.alternateStandardImage)
                if self.imtest(self.spaceTransform):
                    regOptions += "-w %s" % self.spaceTransform
                else:
                    regOptions += "--postmat=%s" % self.spaceTransform

                if self.imtest(self.outputDir + "/native_space/perfusion_calib"):
                    self.runcmd("%s %s -i %s/native_space/perfusion_calib -a %s/native_space/perfusion_calib_standard" %
                              (self.fslProg("applywarp"), regOptions, self.outputDir, self.outputDir))
                else:
                    self.runcmd("%s %s -i %s/native_space/perfusion -a %s/native_space/perfusion_standard" %
                              (self.fslProg("applywarp"), regOptions, self.outputDir, self.outputDir))
                self.runcmd("%s %s -i %s/native_space/arrival -a %s/native_space/arrival_standard" %
                      (self.fslProg("applywarp"), regOptions, self.outputDir, self.outputDir))
        except RuntimeError, e:
            self.warn(str(e), "Options are not valid")

    def __init__(self):
        #self.addListener('inputImage', 'setOutputImage', self.setOutputImage)
        pass


optLabels = {
    'inputImage'           : 'Input image',
    'inversionTimes'       : 'Inversion times',
    'bolusDuration'        : 'Bolus duration',
    'labelling'            : 'Labelling',
    'tagControlPairs'      : 'Data is tag-control pairs',
    'controlFirst'         : 'Data has control first',
    'dataOrder'            : 'Data order (grouped by)',
    'staticTissue'         : 'Static tissue',
    'useStructuralImage'   : 'Structural image',
    'structuralImage'      : '',
    'runBet'               : 'Run BET on Structural image',

    'useSpaceTransform'        : 'Structural to standard space transform',
    'spaceTransform'           : '',
    'useAlternateStandardImage': 'Alternate standard brain image',
    'alternateStandardImage'   : '',

    'outputDir'            : 'Output directory',
    'brainMask'            : 'Optional Brain mask',
    'paramVariance'        : 'Output parameter variance',
    'bolusArrivalTime'     : 'Bolus arrival time',
    't1'                   : 'T1',
    't1b'                  : 'T1b',
    'inversionEffic'       : 'Inversion efficiencey',
    'spatialSmoothing'     : 'Use adaptive spatial smoothing on CBF',
    'inferT1'              : 'Incorporate T1 value uncertainty',
    'includeVascular'      : 'Include macro vascular component',
    'fixBolusDuration'     : 'Fix bolus duration',

    'doCalibration'            : 'Perform calibration',
    'calibrateMode'        : 'Mode',
    'm0Image'              : 'M0 calibration image',
    'useCoilSensitivityRefImage'      : 'Use coil sensitivity reference image',
    'coilSensitivityRefImage'      : '',
    'calibrateGain'        : 'Calibration gain',
    'refTissueType'        : 'Reference tissue type',
    'useRefTissueMask'     : 'Use reference tissue mask',
    'refTissueMask'        : 'Reference tissue mask image',
    'refT1'                : 'Reference T1(s)',
    'refT2'                : 'Reference T2(s)',
    'bloodT2'              : 'Blood T2(s)',
    'seqTR'                : 'Sequence TR(s)',
    'seqTE'                : 'Sequence TE(s)',
}

optTooltips = {
}

aslView = props.VGroup(
    label="Oxford ASL",
    children=(
        props.NotebookGroup((
            props.VGroup(
                label="Data",
                children=(
                    'inputImage',
                    'inversionTimes',
                    'bolusDuration',
                    'labelling',
                    'tagControlPairs',
                    'controlFirst',
                    'dataOrder',
                    conditional('staticTissue','tagControlPairs'),
                    'useStructuralImage',
                    conditional('structuralImage','useStructuralImage'),
                    conditional('runBet', 'useStructuralImage'),
                )
            ),
           props.VGroup(
                label='Analysis',
                children=(
                    'outputDir',
                    'brainMask',
                    'paramVariance',
                    'bolusArrivalTime',
                    't1',
                    't1b',
                    'inversionEffic',
                    'spatialSmoothing',
                    'inferT1',
                    'includeVascular',
                    'fixBolusDuration',
                )
            ),
            props.VGroup(
                label="Registration",
                children=(
                    conditional('useSpaceTransform','useStructuralImage'),
                    conditional('spaceTransform','useSpaceTransform'),
                    conditional('useAlternateStandardImage','useStructuralImage'),
                    conditional('alternateStandardImage','useAlternateStandardImage'),
                )
            ),
            props.VGroup(
                label='Calibration',
                children=(
                    'doCalibration',
                    props.VGroup(
                        children=(
                            'calibrateMode',
                            'm0Image',
                            'useCoilSensitivityRefImage',
                            conditional('coilSensitivityRefImage','useCoilSensitivityRefImage'),
                            'calibrateGain',
                            'refTissueType',
                            't1b',
                            'refT1',
                            'refT2',
                            'bloodT2',
                            'seqTR',
                            'seqTE',
                        ),
                        visibleWhen=lambda i: getattr(i, 'doCalibration')
                    ),
                )
            )
        )),
        props.Button(text='Run', callback=lambda i, b: i.runAsl()),
    )
)

if __name__ == '__main__':

    opts = AslOptions()
    app  = wx.App()
    props.initGUI()
    dlg  = props.buildDialog(None, opts, aslView, optLabels, optTooltips)

    dlg.Show()
    app.MainLoop()
