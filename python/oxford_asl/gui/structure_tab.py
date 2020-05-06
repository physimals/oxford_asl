"""
BASIL GUI for Oxford ASL - Structural data options

Copyright (C) 2020 University of Oxford
"""
import wx
import wx.grid

from .widgets import TabPage

class StructureTab(TabPage):
    """
    Tab page containing options for structural space transformation
    """

    EXISTING_FSLANAT = 0
    NEW_FSLANAT = 1
    INDEP_STRUC = 2
    NONE = 3

    TRANS_MATRIX = 0
    TRANS_IMAGE = 1
    TRANS_FSLANAT = 2

    def __init__(self, app, parent, idx, n):
        TabPage.__init__(self, app, parent, "Structure", idx, n, name="structure")
        self._outdir = "."

        self.section("Structure")

        self.struc_ch = wx.Choice(self, choices=["Existing FSL_ANAT output", "Run FSL_ANAT on structural image", "Independent structural data", "None"])
        self.struc_ch.SetSelection(self.NONE)
        self.struc_ch.Bind(wx.EVT_CHOICE, self.state_changed)
        self.struc_ch.span = 2
        self.pack("Structural data from", self.struc_ch)

        self.fsl_anat_picker = self.file_picker("Existing FSL_ANAT directory", pick_dir=True)
        self.struc_image_picker = self.file_picker("Structural Image")
        self.brain_image_picker = self.file_picker("Brain image", optional=True)

        self.section("Registration")

        self.transform_choices = ["Use matrix", "Use warp image", "Use FSL_ANAT"]

        self.transform_cb = wx.CheckBox(self, label="Transform to standard space")
        self.transform_cb.Bind(wx.EVT_CHECKBOX, self.state_changed)
        self.transform_ch = wx.Choice(self, choices=self.transform_choices)
        self.transform_ch.SetSelection(self.TRANS_FSLANAT)
        self.transform_ch.Bind(wx.EVT_CHOICE, self.state_changed)
        self.transform_picker = wx.FilePickerCtrl(self)
        self.transform_picker.Bind(wx.EVT_FILEPICKER_CHANGED, self.state_changed)
        self.pack("", self.transform_cb, self.transform_ch, self.transform_picker, enable=False)

        self.sizer.AddGrowableCol(2, 1)
        self.SetSizer(self.sizer)
        self.add_next_prev_btn()
        self._update_ui()

    def options(self):
        struc_method = self.struc_ch.GetSelection()
        opts = {
            "_struc_method" : struc_method,
            "s"       : self.struc_image_picker.GetPath() if struc_method in (self.NEW_FSLANAT, self.INDEP_STRUC) else None,
            "sbrain"  : self.brain_image_picker.GetPath() if struc_method == self.INDEP_STRUC and self.brain_image_picker.checkbox.IsChecked() else None,
        }

        if struc_method == self.EXISTING_FSLANAT:
            opts["fslanat"] = self.fsl_anat_picker.GetPath()
        elif struc_method == self.NEW_FSLANAT:
            opts["fslanat"] = "%s/struc.anat" % self._outdir

        # Standard space transform
        if self._have_transform():
            if self._transform_type() == self.TRANS_MATRIX:
                opts["t"] = self.transform_picker.GetPath()
            elif self._transform_type() == self.TRANS_IMAGE:
                opts["warp"] = self.transform_picker.GetPath()
            else:
                # This implies that FSLANAT output is being used, and hence
                # --fslanat is already specified
                pass

        return opts

    def check_options(self, options):
        self._check_image("Structural image", options["s"])
        self._check_image("Brain extracted structural image", options["sbrain"])
        if "t" in options:
            self._check_exists("Standard space transformation matrix", options["t"], can_be_none=False)
        if "warp" in options:
            self._check_image("Warp image", options["warp"], can_be_none=False)
        if "fslanat" in options and self.struc_ch.GetSelection() == self.EXISTING_FSLANAT:
            self._check_exists("FSL_ANAT output", options["fslanat"], can_be_none=False)

    def option_changed(self, options, key, value):
        if key == "i":
            self.set_picker_dir(self.struc_image_picker, value)
            self.set_picker_dir(self.brain_image_picker, value)
            self.set_picker_dir(self.fsl_anat_picker, value)
            self.set_picker_dir(self.transform_picker, value)

        if key == "o":
            self._outdir = value

    def state_changed(self, event=None):
        self._update_ui()
        TabPage.state_changed(self, event)

    def _update_ui(self):
        mode = self.struc_ch.GetSelection()
        self.fsl_anat_picker.Enable(mode == self.EXISTING_FSLANAT)
        self.struc_image_picker.Enable(mode in (self.NEW_FSLANAT, self.INDEP_STRUC))

        self.brain_image_picker.checkbox.Enable(mode == self.INDEP_STRUC)
        self.brain_image_picker.Enable(mode == self.INDEP_STRUC and self.brain_image_picker.checkbox.IsChecked())

        # Only offer FSL_ANAT transform option if we are using FSL_ANAT
        sel = self.transform_ch.GetSelection()
        if mode == self.INDEP_STRUC:
            if sel == self.TRANS_FSLANAT:
                sel = self.TRANS_MATRIX
            choices = 2
        else:
            choices = 3
        self.transform_ch.Enable(False)
        self.transform_ch.Clear()
        self.transform_ch.AppendItems(self.transform_choices[:choices])
        self.transform_ch.SetSelection(sel)
        self.transform_ch.Enable(self._have_transform())

        self.transform_picker.Enable(self._have_transform() and self._transform_type() != self.TRANS_FSLANAT)

    def _have_transform(self):
        return self.transform_cb.IsChecked()

    def _transform_type(self):
        return self.transform_ch.GetSelection()
