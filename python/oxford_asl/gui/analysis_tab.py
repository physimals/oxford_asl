"""
BASIL GUI for Oxford ASL - Analysis / model fitting options

Copyright (C) 2020 University of Oxford
"""
import os

import wx

from . import OptionError
from .widgets import TabPage, WhitePaperCompatibility

class AnalysisTab(TabPage):
    """
    Tab page containing data analysis options
    """

    def __init__(self, app, parent, idx, n):
        TabPage.__init__(self, app, parent, "Analysis", idx, n)

        self.distcorr_choices = ["Fieldmap", "Calibration image"]

        self.section("Basic analysis options")

        self.outdir_picker = self.file_picker("Output Directory", pick_dir=True)
        self.mask_picker = self.file_picker("User-specified brain Mask", optional=True)

        self.section("White paper compatibility")
        #self.wp_reset = wx.Button(self, label="Reset to white paper defaults")
       
        #self.wp_cb = self.checkbox("Analysis which conforms to 'White Paper' (Alsop et al 2014)")

        self.section("Initial parameter values")

        self.bat_num = self.number("Arterial Transit Time (s)", minval=0, maxval=2.5, initial=1.3)
        self.t1_num = self.number("T1 (s)", minval=0, maxval=3, initial=1.3)
        self.t1b_num = self.number("T1b (s)", minval=0, maxval=3, initial=1.65)
        self.ie_num = self.number("Inversion Efficiency", minval=0, maxval=1, initial=0.85)

        self.section("Analysis Options")

        self.spatial_cb = self.checkbox("Adaptive spatial regularization on perfusion", initial=True)
        self.infer_t1_cb = self.checkbox("Incorporate T1 value uncertainty")
        self.macro_cb = self.checkbox("Include macro vascular component")
        self.fixbolus_cb = self.checkbox("Fix label duration", initial=True)

        self.pv_cb = self.checkbox("Partial Volume Correction")
        self.mc_cb = self.checkbox("Motion Correction")

        self.sizer.AddGrowableCol(1, 1)
        self.SetSizer(self.sizer)
        self.add_next_prev_btn()

    def options(self):
        opts = {
            "o"           : self.outdir_picker.GetPath(),
            "m"           : self.mask_picker.GetPath() if self.mask_picker.checkbox.IsChecked() else None,
            #"wp"          : self.wp_cb.IsChecked(),
            "bat"         : self.bat_num.GetValue(),
            "t1"          : self.t1_num.GetValue(),
            "t1b"         : self.t1b_num.GetValue(),
            "alpha"       : self.ie_num.GetValue(),
            "spatial"     : self.spatial_cb.IsChecked(),
            "fixbolus"    : self.fixbolus_cb.IsChecked(),
            "mc"          : self.mc_cb.IsChecked(),
            "infert1"     : self.infer_t1_cb.IsChecked(),
            "pvcorr"      : self.pv_cb.IsChecked(),
            "artoff"      : not self.macro_cb.IsChecked()
        }

        #if opts["wp"]:
        #    opts.pop("t1")
        #    opts.pop("bat")
        return opts

    def check_options(self, options):
        self._check_image("Analysis mask", options["m"])

        outdir = options["o"]
        if not outdir:
            raise OptionError("Output directory not specified")
        elif os.path.exists(outdir) and not os.path.isdir(outdir):
            raise OptionError("Output directory already exists and is a file")

    def option_changed(self, options, key, value):
        if key == "i":
            self.set_picker_dir(self.outdir_picker, value)
            self.set_picker_dir(self.mask_picker, value)

        if key == "wp":
            if value:
                self.t1_num.SetValue(1.65)
                self.bat_num.SetValue(0)
            else:
                self.t1_num.SetValue(1.3)
                self.bat_num.SetValue(1.3)

        if key == "casl":
            if value:
                self.bat_num.SetValue(1.3)
                self.ie_num.SetValue(0.85)
            else:
                self.bat_num.SetValue(0.7)
                self.ie_num.SetValue(0.98)
