import wx
import wx.grid

from .widgets import TabPage

class AslCalibration(TabPage):
    """ 
    Tab page containing calibration options
    """

    def __init__(self, parent, idx, n):
        TabPage.__init__(self, parent, "Calibration", idx, n)

        self.calib_cb = self.checkbox("Enable Calibration", bold=True, handler=self.calib_changed)

        self.calib_image_picker = self.file_picker("Calibration Image")
        self.m0_type_ch = self.choice("M0 Type", choices=["Proton Density (long TR)", "Saturation Recovery"])

        self.seq_tr_num = self.number("Sequence TR (s)", min=0,max=10,initial=6)
        self.calib_gain_num = self.number("Calibration Gain", min=0,max=5,initial=1)
        self.calib_mode_ch = self.choice("Calibration mode", choices=["Reference Region", "Voxelwise"])

        self.section("Reference tissue")

        self.ref_tissue_type_ch = self.choice("Type", choices=["CSF", "WM", "GM", "None"], handler=self.ref_tissue_type_changed)
        self.ref_tissue_mask_picker = self.file_picker("Mask", optional=True)
        self.ref_t1_num = self.number("Reference T1 (s)", min=0,max=5,initial=4.3)
        self.seq_te_num = self.number("Sequence TE (ms)", min=0,max=30,initial=0)
        self.ref_t2_num = self.number("Reference T2 (ms)", min=0,max=1000,initial=750, step=10)
        self.blood_t2_num = self.number("Blood T2 (ms)", min=0,max=1000,initial=150, step=10)
        self.coil_image_picker = self.file_picker("Coil Sensitivity Image", optional=True)

        self.sizer.AddGrowableCol(2, 1)
        self.SetSizer(self.sizer)
        self.next_prev()

    def calib(self): return self.calib_cb.IsChecked()
    def m0_type(self): return self.m0_type_ch.GetSelection()
    def seq_tr(self): return self.seq_tr_num.GetValue()
    def seq_te(self): return self.seq_te_num.GetValue()
    def calib_image(self): return self.calib_image_picker.GetPath()
    def calib_gain(self): return self.calib_gain_num.GetValue()
    def calib_mode(self): return self.calib_mode_ch.GetSelection()
    def ref_tissue_type(self): return self.ref_tissue_type_ch.GetSelection()
    def ref_tissue_type_name(self): return self.ref_tissue_type_ch.GetString(self.ref_tissue_type())
    def ref_tissue_mask(self): 
        if self.ref_tissue_mask_picker.checkbox.IsChecked():
            return self.ref_tissue_mask_picker.GetPath()
        else:
            return None
    def ref_t1(self): return self.ref_t1_num.GetValue()
    def ref_t2(self): return self.ref_t2_num.GetValue()
    def blood_t2(self): return self.blood_t2_num.GetValue()
    def coil_image(self): 
        if self.coil_image_picker.checkbox.IsChecked(): return self.coil_image_picker.GetPath()
        else: return None

    def ref_tissue_type_changed(self, event):
        if self.ref_tissue_type() == 0: # CSF
            self.ref_t1_num.SetValue(4.3)
            self.ref_t2_num.SetValue(750)
        elif self.ref_tissue_type() == 1: # WM
            self.ref_t1_num.SetValue(1.0)
            self.ref_t2_num.SetValue(50)
        elif self.ref_tissue_type() == 2: # GM
            self.ref_t1_num.SetValue(1.3)
            self.ref_t2_num.SetValue(100)
        self.update()

    def calib_changed(self, event):
        self.distcorr.calib_changed(self.calib())
        self.update()

    def wp_changed(self, wp):
        self.update()

    def update(self, event=None):
        enable = self.calib()
        self.m0_type_ch.Enable(enable)
        self.seq_tr_num.Enable(enable and self.m0_type() == 0)
        self.calib_image_picker.Enable(enable)
        self.calib_gain_num.Enable(enable)
        self.coil_image_picker.checkbox.Enable(enable)
        if self.analysis.wp(): self.calib_mode_ch.SetSelection(1)
        self.calib_mode_ch.Enable(enable and not self.analysis.wp())
        self.ref_tissue_type_ch.Enable(enable and self.calib_mode() == 0)

        if self.ref_tissue_type() == 3:
            # Ref tissue = None - enforce mask
            self.ref_tissue_mask_picker.checkbox.Enable(False)
            self.ref_tissue_mask_picker.checkbox.SetValue(enable and self.calib_mode() == 0)
            self.ref_tissue_mask_picker.Enable(enable and self.calib_mode() == 0)
        else:
            self.ref_tissue_mask_picker.checkbox.Enable(enable and self.calib_mode() == 0)
        self.ref_tissue_mask_picker.Enable(enable and self.ref_tissue_mask_picker.checkbox.IsChecked() and self.calib_mode() == 0)
        
        self.coil_image_picker.checkbox.Enable(enable and self.calib_mode() == 0)
        self.coil_image_picker.Enable(enable and self.calib_mode() == 0 and self.coil_image_picker.checkbox.IsChecked())
        self.seq_te_num.Enable(enable and self.calib_mode() == 0)
        self.blood_t2_num.Enable(enable and self.calib_mode() == 0)
        self.ref_t1_num.Enable(enable and self.calib_mode() == 0)
        self.ref_t2_num.Enable(enable and self.calib_mode() == 0)
        TabPage.update(self)
