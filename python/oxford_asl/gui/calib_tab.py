
"""
BASIL GUI for Oxford ASL - Calibration options

Copyright (C) 2020 University of Oxford
"""
from .widgets import TabPage

# T1 and T2 defaults for different reference tissue types
REF_TISSUE_DEFAULTS = {
    "csf" : (4.3, 750),
    "wm" : (1.0, 50),
    "gm" : (1.3, 100),
}

class CalibTab(TabPage):
    """
    Tab page containing calibration options
    """
    def __init__(self, app, parent, idx, n):
        TabPage.__init__(self, app, parent, "Calibration", idx, n)

        self.calib_cb = self.checkbox("Enable Calibration", bold=True)

        self.calib_image_picker = self.file_picker("Calibration Image")
        self.m0_type_ch = self.choice("M0 Type", choices=["Proton Density (long TR)", "Saturation Recovery"])

        self.seq_tr_num = self.number("Sequence TR (s)", minval=0, maxval=10, initial=6)
        self.calib_gain_num = self.number("Calibration Gain", minval=0, maxval=5, initial=1)
        self.calib_mode_ch = self.choice("Calibration mode", choices=["Reference Region", "Voxelwise"])

        self.section("Reference tissue")
        self.ref_tissue_type_ch = self.choice("Type", choices=["CSF", "WM", "GM", "None"])
        self.ref_tissue_mask_picker = self.file_picker("Mask", optional=True)
        self.ref_t1_num = self.number("Reference T1 (s)", minval=0, maxval=5, initial=4.3)
        self.seq_te_num = self.number("Sequence TE (ms)", minval=0, maxval=30, initial=0)
        self.ref_t2_num = self.number("Reference T2 (ms)", minval=0, maxval=1000, initial=750, step=10)
        self.blood_t2_num = self.number("Blood T2 (ms)", minval=0, maxval=1000, initial=150, step=10)
        self.coil_image_picker = self.file_picker("Reference Image for sensitivity correction", optional=True)

        self.sizer.AddGrowableCol(2, 1)
        self.SetSizer(self.sizer)
        self.add_next_prev_btn()
        self._update_ui()

    def options(self):
        # FIXME saturation recovery
        opts = {
            "_calib_enabled" : self.calib_cb.IsChecked(),
        }
        if opts["_calib_enabled"]:
            opts.update({
                "c"       : self.calib_image_picker.GetPath(),
                "cmethod" : "single" if self.calib_mode_ch.GetSelection() == 0 else "voxel",
                "tr"      : self.seq_tr_num.GetValue(),
                "cgain"   : self.calib_gain_num.GetValue(),
                "cref"    : self.coil_image_picker.GetPath() if self.coil_image_picker.checkbox.IsChecked() else None,
            })

            if opts["cmethod"] == "single":
                opts.update({
                    "tissref" : self.ref_tissue_type_ch.GetString(self.ref_tissue_type_ch.GetSelection()).lower(),
                    "csf"     : self.ref_tissue_mask_picker.GetPath() if self.ref_tissue_mask_picker.checkbox.IsChecked() else None,
                    "t1csf"   : self.ref_t1_num.GetValue(),
                    "t2csf"   : self.ref_t2_num.GetValue(),
                    "t2bl"    : self.blood_t2_num.GetValue(),
                    "te"      : self.seq_te_num.GetValue(),
                })

        return opts

    def check_options(self, options):
        if "c" in options:
            self._check_image("Calibration image", options["c"], can_be_none=False)
            self._check_image("Coil sensitivity reference image", options["cref"])
            self._check_image("Calibration reference tissue mask", options.get("csf", None))

    def state_changed(self, event=None):
        self._update_ui()
        TabPage.state_changed(self, event)

    def _update_ui(self):
        calib_enabled = self.calib_cb.IsChecked()
        refregion = calib_enabled and self.calib_mode_ch.GetSelection() == 0

        self.m0_type_ch.Enable(calib_enabled)
        self.seq_tr_num.Enable(calib_enabled and self.m0_type_ch.GetSelection() == 0)
        self.calib_image_picker.Enable(calib_enabled)
        self.calib_gain_num.Enable(calib_enabled)
        self.coil_image_picker.checkbox.Enable(calib_enabled)
        self.calib_mode_ch.Enable(calib_enabled)

        self.ref_tissue_type_ch.Enable(refregion)
        if refregion and self.ref_tissue_type_ch.GetSelection() == 3:
            # Ref tissue = None - must have user-supplied mask
            self.ref_tissue_mask_picker.checkbox.Enable(False)
            self.ref_tissue_mask_picker.checkbox.SetValue(True)

        self.ref_tissue_mask_picker.checkbox.Enable(refregion)
        self.ref_tissue_mask_picker.Enable(refregion and self.ref_tissue_mask_picker.checkbox.IsChecked())

        self.coil_image_picker.checkbox.Enable(refregion)
        self.coil_image_picker.Enable(refregion and self.coil_image_picker.checkbox.IsChecked())
        self.seq_te_num.Enable(refregion)
        self.blood_t2_num.Enable(refregion)
        self.ref_t1_num.Enable(refregion)
        self.ref_t2_num.Enable(refregion)

    def option_changed(self, options, key, value):
        if key == "i":
            self.set_picker_dir(self.calib_image_picker, value)
            self.set_picker_dir(self.coil_image_picker, value)
            self.set_picker_dir(self.ref_tissue_mask_picker, value)

        if key == "tissref":
            ref_t1, ref_t2 = REF_TISSUE_DEFAULTS.get(value, (None, None))
            if ref_t1 is not None:
                self.ref_t1_num.SetValue(ref_t1)
                self.ref_t2_num.SetValue(ref_t2)

        if key == "wp" and value:
            # Set wp defaults
            self.calib_cb.SetValue(True)
            self.calib_mode_ch.SetSelection(1)
            self._update_ui()
