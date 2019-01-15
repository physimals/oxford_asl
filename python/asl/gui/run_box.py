import sys
import os
import traceback
import tempfile
import shutil
import shlex
import subprocess
from threading import Thread

import nibabel as nib

import wx
from wx.lib.pubsub import pub

class OptionError(RuntimeError):
    pass

class Mkdir:
    def __init__(self, dirname):
        self.dirname = dirname

    def run(self):
        if not os.path.exists(self.dirname):
            os.makedirs(self.dirname)
        return 0

class FslCmd:
    def __init__(self, cmd):
        script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
        fsldevdir = os.path.join(os.environ.get("FSLDEVDIR", ""), "bin")
        fsldir = os.path.join(os.environ.get("FSLDIR", ""), "bin")

        self.cmd = cmd
        for d in (script_dir, fsldevdir, fsldir):
            if os.path.exists(os.path.join(d, cmd)):
                self.cmd = os.path.join(d, cmd)
                break
            
    def add(self, opt, val=None):
        if val is not None:
            self.cmd += " %s=%s" % (opt, str(val))
        else:
            self.cmd += " %s" % opt

    def write_output(self, line):
        wx.CallAfter(pub.sendMessage, "run_stdout", line=line)

    def run(self):
        self.write_output(self.cmd + "\n")
        args = shlex.split(self.cmd)
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        while 1:
            retcode = p.poll() #returns None while subprocess is running
            line = p.stdout.readline().decode('utf-8')
            self.write_output(line)
            if retcode is not None: break
        self.write_output("\nReturn code: %i\n\n" % retcode)
        return retcode

    def __str__(self): return self.cmd

class CmdRunner(Thread):
    def __init__(self, cmds, done_cb):
        Thread.__init__(self)
        self.cmds = cmds
        self.done_cb = done_cb

    def run(self):
        ret = -1
        try:
            for cmd in self.cmds:
                ret = -1
                ret = cmd.run()
                if ret != 0:
                    break
        finally:
            wx.CallAfter(pub.sendMessage, "run_finished", retcode=ret)

class AslRun(wx.Frame):
    """
    Determines the commands to run and displays them in a window
    """

    # The options we need to pass to oxford_asl for various data orderings
    order_opts = {"trp" : "--ibf=rpt --iaf=diff", 
                  "trp,tc" : "--ibf=rpt --iaf=tcb", 
                  "trp,ct" : "--ibf=rpt --iaf=ctb",
                  "rtp" : "--ibf=tis --iaf=diff",
                  "rtp,tc" : "--ibf=tis --iaf=tcb",
                  "rtp,ct" : "--ibf=tis --iaf=ctb",
                  "ptr,tc" : "--ibf=rpt --iaf=tc",
                  "ptr,ct" : "--ibf=rpt --iaf=ct",
                  "prt,tc" : "--ibf=tis --iaf=tc",
                  "prt,ct" : "--ibf=tis --iaf=ct"}

    def __init__(self, parent, run_btn, run_label):
        wx.Frame.__init__(self, parent, title="Run", size=(600, 400), style=wx.DEFAULT_FRAME_STYLE)

        self.run_seq = None
        self.run_btn = run_btn
        self.run_btn.Bind(wx.EVT_BUTTON, self.dorun)
        self.run_label = run_label
        self.preview_data = None
    
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.output_text = wx.TextCtrl(self, style=wx.TE_READONLY | wx.TE_MULTILINE)
        font = wx.Font(8, wx.TELETYPE, wx.NORMAL, wx.NORMAL)
        self.output_text.SetFont(font)
        self.sizer.Add(self.output_text, 1, flag=wx.EXPAND)
            
        self.SetSizer(self.sizer)
        self.Bind(wx.EVT_CLOSE, self.close)
        pub.subscribe(self.write_output, "run_stdout")
        pub.subscribe(self.finished, "run_finished")

    def write_output(self, line):
        self.output_text.AppendText(line)

    def close(self, _):
        self.Hide()

    def finished(self, retcode):
        if retcode != 0:
            self.write_output("\nWARNING: command failed\n")
        self.update()

    def dorun(self, _):
        if self.run_seq: 
            self.Show()
            self.Raise()
            self.output_text.Clear()
            self.run_btn.Enable(False)
            self.run_label.SetForegroundColour(wx.Colour(0, 0, 128))
            self.run_label.SetLabel("Running - Please Wait")
            runner = CmdRunner(self.run_seq, self.finished)
            runner.start()

    def update(self):
        """
        Get the sequence of commands and enable the run button if options are valid. Otherwise
        display the first error in the status label
        """
        self.run_seq = None
        try:
            self.run_seq = self.get_run_sequence()
            self.run_label.SetForegroundColour(wx.Colour(0, 128, 0))
            self.run_label.SetLabel("Ready to Go")
            self.run_btn.Enable(True)
        except (OptionError, nib.filebasedimages.ImageFileError) as e:
            self.run_btn.Enable(False)
            self.run_label.SetForegroundColour(wx.Colour(255, 0, 0))
            self.run_label.SetLabel(str(e))
        except:
            # Any other exception is a program bug - report it to STDERR
            self.run_btn.Enable(False)
            self.run_label.SetForegroundColour(wx.Colour(255, 0, 0))
            self.run_label.SetLabel("Unexpected error - see console and report as a bug")
            traceback.print_exc(sys.exc_info()[1])

    def check_exists(self, label, fname):
        if not os.path.exists(fname):
            raise OptionError("%s - no such file or directory" % label)

    def get_preview_data(self):
        """
        Run ASL_FILE for perfusion weighted image - just for the preview
        """
        infile = self.input.data()
        if infile == "":
            # Don't bother if we have not input file yet!
            return None
        
        tempdir = tempfile.mkdtemp()
        self.preview_data = None
        try:
            meanfile = "%s/mean.nii.gz" % tempdir
            cmd = FslCmd("asl_file")
            cmd.add('--data="%s"' % self.input.data())
            cmd.add("--ntis=%i" % self.input.ntis())
            cmd.add('--mean="%s"' % meanfile)
            cmd.add(" ".join(self.get_data_order_options()))
            cmd.run()
            img = nib.load(meanfile)
            return img.get_data()
        except:
            traceback.print_exc()
            return None
        finally:
            shutil.rmtree(tempdir)

    def get_data_order_options(self):
        """
        Check data order is supported and return the relevant options
        """
        order, tagfirst = self.input.data_order()
        diff_opt = ""
        if self.input.tc_pairs(): 
            if tagfirst: order += ",tc"
            else: order += ",ct"
            diff_opt = "--diff"
        if order not in self.order_opts:
            raise OptionError("This data ordering is not supported by ASL_FILE")
        else: 
            return self.order_opts[order], diff_opt

    def get_run_sequence(self):
        """
        Get the sequence of commands for the selected options, throwing exception
        if any problems are found (e.g. files don't exist, mandatory options not specified)

        Exception text is reported by the GUI
        """
        run = []

        # Check input file exists, is an image and the TIs/repeats/TC pairs is consistent
        self.check_exists("Input data", self.input.data())
        img = nib.load(self.input.data())
        if len(img.shape) != 4:
            raise OptionError("Input data is not a 4D image")
        nvols = img.shape[3]

        N = self.input.ntis()
        if self.input.tc_pairs(): N *= 2
        if nvols % N != 0:
            self.input.nrepeats_label.SetLabel("<Invalid>")
            raise OptionError("Input data contains %i volumes - not consistent with %i TIs and TC pairs=%s" % (img.shape[3], self.input.ntis(), self.input.tc_pairs()))
        else:
            self.input.nrepeats_label.SetLabel("%i" % (nvols / N))
            self.preview.order_preview.n_tis = self.input.ntis()
            self.preview.order_preview.n_repeats = nvols / N
            self.preview.order_preview.tc_pairs = self.input.tc_pairs()
            self.preview.order_preview.tagfirst = self.input.tc_ch.GetSelection() == 0
            self.preview.order_preview.Refresh()

        # Build OXFORD_ASL command 
        outdir = self.analysis.outdir()
        if os.path.exists(outdir) and not os.path.isdir(outdir):
            raise OptionError("Output directory already exists and is a file")
        run.append(Mkdir(outdir))

        # Input data
        cmd = FslCmd("oxford_asl")
        cmd.add(' -i "%s"' % self.input.data())
        cmd.add(self.get_data_order_options()[0])
        cmd.add("--tis %s" % ",".join(["%.2f" % v for v in self.input.tis()]))
        cmd.add("--bolus %s" % ",".join(["%.2f" % v for v in self.input.bolus_dur()]))
        if self.input.labelling() == 1: 
            cmd.add("--casl")
        if self.input.readout() == 1:
            # 2D multi-slice readout - must give dt in seconds
            cmd.add("--slicedt %.5f" % (self.input.time_per_slice() / 1000))
            if self.input.multiband():
                cmd.add("--sliceband %i" % self.input.slices_per_band())

        # Structure - may require FSL_ANAT to be run
        fsl_anat_dir = self.structure.existing_fsl_anat()
        struc_image = self.structure.struc_image()
        if fsl_anat_dir is not None:
            # Have an existing FSL_ANAT directory
            self.check_exists("FSL_ANAT", fsl_anat_dir)
            cmd.add('--fslanat="%s"' % fsl_anat_dir)
        elif self.structure.run_fsl_anat():
            # FIXME set this up and pass in the dir using --fslanat
            self.check_exists("Structural image", struc_image)
            fsl_anat = FslCmd("fsl_anat")
            fsl_anat.add('-i "%s"' % struc_image)
            fsl_anat.add('-o "%s/struc"' % outdir)
            run.append(fsl_anat)
            cmd.add('--fslanat="%s/struc.anat"' % outdir)
        elif struc_image is not None:
            # Providing independent structural data
            self.check_exists("Structural image", struc_image)
            cp = FslCmd("imcp")
            cp.add('"%s"' % struc_image)
            cp.add('"%s/structural_head"' % outdir)
            run.append(cp)
            cmd.add('-s "%s/structural_head"' % outdir)

            # Brain image can be provided or can use BET
            brain_image = self.structure.struc_image_brain()
            if brain_image is not None:
                self.check_exists("Structural brain image", brain_image)
                cp = FslCmd("imcp")
                cp.add('"%s"' % brain_image)
                cp.add('"%s/structural_brain"' % outdir)
                run.append(cp)
            else:
                bet = FslCmd("bet")
                bet.add('"%s"' % struc_image)
                bet.add('"%s/structural_brain"' % outdir)
                run.append(bet)
            cmd.add('--sbrain "%s/structural_brain"' % outdir)
        else:
            # No structural data
            pass
        # Structure transform
        if self.structure.transform():
            if self.structure.transform_type() == self.structure.TRANS_MATRIX:
                self.check_exists("Transformation matrix", self.structure.transform_file())
                cmd.add('-t "%s"' % self.structure.transform_file())
            elif self.structure.transform_type() == self.structure.TRANS_IMAGE:
                self.check_exists("Warp image", self.structure.transform_file())
                cmd.add('--warp "%s"' % self.structure.transform_file())
            else:
                # This implies that FSLANAT output is being used, and hence
                # --fslanat is already specified
                pass 

        # Calibration - do this via oxford_asl rather than calling asl_calib separately
        if self.calibration.calib():
            self.check_exists("Calibration image", self.calibration.calib_image())
            cmd.add('-c "%s"' % self.calibration.calib_image())
            if self.calibration.m0_type() == 0:
                #calib.add("--mode longtr")
                cmd.add("--tr %.2f" % self.calibration.seq_tr())
            else:
                raise OptionError("Saturation recovery not supported by oxford_asl")
                #calib.add("--mode satrevoc")
                #calib.add("--tis %s" % ",".join([str(v) for v in self.input.tis()]))
                # FIXME change -c option in sat recov mode?

            cmd.add("--cgain %.2f" % self.calibration.calib_gain())
            if self.calibration.calib_mode() == 0:
                cmd.add("--cmethod single")
                cmd.add("--tissref %s" % self.calibration.ref_tissue_type_name().lower())
                cmd.add("--te %.2f" % self.calibration.seq_te())
                cmd.add("--t1csf %.2f" % self.calibration.ref_t1())
                cmd.add("--t2csf %.2f" % self.calibration.ref_t2())
                cmd.add("--t2bl %.2f" % self.calibration.blood_t2())
                if self.calibration.ref_tissue_mask() is not None:
                    self.check_exists("Calibration reference tissue mask", self.calibration.ref_tissue_mask())
                    cmd.add('--csf "%s"' % self.calibration.ref_tissue_mask())
                if self.calibration.coil_image() is not None:
                    self.check_exists("Coil sensitivity reference image", self.calibration.coil_image())
                    cmd.add('--cref "%s"' % self.calibration.coil_image())
            else:
                cmd.add("--cmethod voxel")
      
        # Distortion correction
        if self.distcorr.distcorr():
            if self.distcorr.distcorr_type() == self.distcorr.FIELDMAP:
                # Fieldmap image
                fmap = self.distcorr.fmap()
                self.check_exists("Fieldmap image", fmap)
                cmd.add('--fmap="%s"' % fmap)
                fmap_mag = self.distcorr.fmap_mag()
                self.check_exists("Fieldmap magnitude image", fmap_mag)
                cmd.add('--fmapmag="%s"' % fmap_mag)
                fmap_be = self.distcorr.fmap_mag_be()
                if fmap_be is not None:
                    self.check_exists("Brain-extracted fieldmap magnitude image", fmap_be)
                    cmd.add('--fmapmagbrain="%s"' % fmap_be)
            else:
                # Calibration image
                calib = self.distcorr.calib()
                self.check_exists("Phase encode reversed calibration image", calib)
                cmd.add('--cblip="%s"' % calib)

            # Generic options - note that echospacing is ms in the GUI but oxford_asl requires it in seconds
            cmd.add("--echospacing=%.5f" % (self.distcorr.echosp() / 1000))
            cmd.add("--pedir=%s" % self.distcorr.pedir())
  
        # Analysis options
        if self.analysis.wp(): 
            cmd.add("--wp")
        else: 
            cmd.add("--t1 %.2f" % self.analysis.t1())
            cmd.add("--bat %.2f" % self.analysis.bat())
        cmd.add("--t1b %.2f" % self.analysis.t1b())
        cmd.add("--alpha %.2f" % self.analysis.ie())
        if self.analysis.fixbolus(): cmd.add("--fixbolus")
        if self.analysis.spatial(): cmd.add("--spatial")
        if self.analysis.mc(): cmd.add("--mc")
        if self.analysis.infer_t1(): cmd.add("--infert1")
        if self.analysis.pv(): cmd.add("--pvcorr")
        if not self.analysis.macro(): cmd.add("--artoff")
        if self.analysis.mask() is not None:
            self.check_exists("Analysis mask", self.analysis.mask())
            cmd.add('-m "%s"' % self.analysis.mask())

        # Output dir
        if outdir == "": 
            raise OptionError("Output directory not specified")
        cmd.add('-o "%s"' % outdir)

        run.append(cmd)
        return run
