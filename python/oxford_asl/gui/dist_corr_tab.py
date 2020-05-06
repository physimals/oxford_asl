"""
BASIL GUI for Oxford ASL - Calibration options

Copyright (C) 2020 University of Oxford
"""

from .widgets import TabPage

class DistCorrTab(TabPage):
    """
    Tab page containing distortion correction options
    """

    # Distortion correction choices
    FIELDMAP, CALIB_IMAGE = 0, 1
    DISTCORR_CHOICES = ["Fieldmap", "Calibration image"]

    def __init__(self, app, parent, idx, n):
        TabPage.__init__(self, app, parent, "Distortion Correction", idx, n, name="distcorr")

        self.section("Distortion Correction")

        # Calibration image options
        self.distcorr_ch = self.choice("Apply distortion correction",
                                       self.DISTCORR_CHOICES[:1],
                                       initial=self.FIELDMAP, optional=True)

        # Calib image options
        self.section("Calibration Image Mode")
        self.cblip_picker = self.file_picker("Phase-encode-reversed calibration image")

        # Fieldmap options
        self.section("Fieldmap Mode")
        self.fmap_picker = self.file_picker("Fieldmap image (in rad/s)")
        self.fmap_mag_picker = self.file_picker("Fieldmap magnitude image")
        self.fmap_be_picker = self.file_picker("Brain-extracted magnitude image", optional=True)

        # General options
        self.section("General")
        self.echosp_num = self.number("Effective EPI echo spacing (ms)", minval=0, maxval=10, digits=5, initial=1.0)
        self.pedir_ch = self.choice("Phase encoding direction", choices=["x", "y", "z", "-x", "-y", "-z"])

        self.sizer.AddGrowableCol(1, 1)
        self.SetSizer(self.sizer)
        self.add_next_prev_btn()

    def options(self):
        if self.distcorr_ch.checkbox.IsChecked():
            opts = {
                "echospacing" : self.echosp_num.GetValue() / 1000,
                "pedir"       : self.pedir_ch.GetStringSelection(),
            }

            if self.distcorr_ch.GetSelection() == self.FIELDMAP:
                opts.update({
                    "fmap" : self.fmap_picker.GetPath(),
                    "fmapmag" : self.fmap_mag_picker.GetPath(),
                    "fmapmagbrain" : self.fmap_be_picker.GetPath() if self.fmap_be_picker.checkbox.IsChecked() else None,
                })

            else:
                opts["cblip"] = self.cblip_picker.GetPath()
            return opts
        else:
            return {}

    def check_options(self, options):
        if "fmap" in options:
            self._check_image("Fieldmap image", options["fmap"], can_be_none=False)
            self._check_image("Fieldmap magnitude image", options["fmapmag"], can_be_none=False)
            self._check_image("Brain-extracted fieldmap magnitude image", options["fmapmagbrain"])
        elif "cblip" in options:
            self._check_image("Phase encode reversed calibration image", options["cblip"], can_be_none=False)

    def state_changed(self, event=None):
        distcorr_enabled = self.distcorr_ch.checkbox.IsChecked()
        self.distcorr_ch.Enable(distcorr_enabled)
        self.pedir_ch.Enable(distcorr_enabled)
        self.echosp_num.Enable(distcorr_enabled)

        fmap = distcorr_enabled and self.distcorr_ch.GetSelection() == self.FIELDMAP
        self.fmap_picker.Enable(fmap)
        self.fmap_mag_picker.Enable(fmap)
        self.fmap_be_picker.Enable(fmap and self.fmap_be_picker.checkbox.IsChecked())
        self.cblip_picker.Enable(distcorr_enabled and not fmap)
        TabPage.state_changed(self, event)

    def option_changed(self, options, key, value):
        if key == "i":
            self.set_picker_dir(self.fmap_picker, value)
            self.set_picker_dir(self.fmap_mag_picker, value)
            self.set_picker_dir(self.fmap_be_picker, value)
            self.set_picker_dir(self.cblip_picker, value)

        # If calibration is enabled, add the phase-reversed calibration image
        # as an option for distortion correction
        if key == "_calib_enabled":
            if value:
                distcorr_choices = self.DISTCORR_CHOICES
                default_choice = 1
            else:
                distcorr_choices = self.DISTCORR_CHOICES[:1]
                default_choice = 0

            self.distcorr_ch.Enable(False)
            self.distcorr_ch.Clear()
            self.distcorr_ch.AppendItems(distcorr_choices)
            self.distcorr_ch.SetSelection(default_choice)
