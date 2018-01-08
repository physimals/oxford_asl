import wx
import wx.grid

from .widgets import TabPage

class AslAnalysis(TabPage):
    """
    Tab page containing data analysis options
    """

    def __init__(self, parent, idx, n):
        TabPage.__init__(self, parent, "Analysis", idx, n)

        self.distcorr_choices = ["Fieldmap", "Calibration image"]

        self.section("Basic analysis options")

        self.outdir_picker = self.file_picker("Output Directory", dir=True)
        self.mask_picker = self.file_picker("Brain Mask", optional=True)
        self.wp_cb = self.checkbox("Analysis which conforms to 'White Paper' (Alsop et al 2014)", handler=self.wp_changed)

        self.section("Initial parameter values")

        self.bat_num = self.number("Arterial Transit Time (s)", min=0,max=2.5,initial=1.3)
        self.t1_num = self.number("T1 (s)", min=0,max=3,initial=1.3)
        self.t1b_num = self.number("T1b (s)", min=0,max=3,initial=1.65)
        self.ie_num = self.number("Inversion Efficiency", min=0,max=1,initial=0.85)
        
        self.section("Analysis Options")

        self.spatial_cb = self.checkbox("Adaptive spatial regularization on perfusion", initial=True)
        self.infer_t1_cb = self.checkbox("Incorporate T1 value uncertainty")
        self.macro_cb = self.checkbox("Include macro vascular component")
        self.fixbolus_cb = self.checkbox("Fix label duration", initial=True)

        self.pv_cb = self.checkbox("Partial Volume Correction")
        self.mc_cb = self.checkbox("Motion Correction")

        self.sizer.AddGrowableCol(1, 1)
        #sizer.AddGrowableRow(5, 1)
        self.SetSizer(self.sizer)
        self.next_prev()

    def outdir(self): return self.outdir_picker.GetPath()
    def mask(self): 
        if self.mask_picker.checkbox.IsChecked(): return self.mask_picker.GetPath()
        else: return None
    def wp(self): return self.wp_cb.IsChecked()
    def bat(self): return self.bat_num.GetValue()
    def t1(self): return self.t1_num.GetValue()
    def t1b(self): return self.t1b_num.GetValue()
    def ie(self): return self.ie_num.GetValue()
    def spatial(self): return self.spatial_cb.IsChecked()
    def infer_t1(self): return self.infer_t1_cb.IsChecked()
    def macro(self): return self.macro_cb.IsChecked()
    def fixbolus(self): return self.fixbolus_cb.IsChecked()
    def pv(self): return self.pv_cb.IsChecked()
    def mc(self): return self.mc_cb.IsChecked()

    def update(self, event=None):
        self.mask_picker.Enable(self.mask_picker.checkbox.IsChecked())
        self.t1_num.Enable(not self.wp())
        self.bat_num.Enable(not self.wp())
        TabPage.update(self)

    def wp_changed(self, event):
        if self.wp():
            self.t1_num.SetValue(1.65)
            self.bat_num.SetValue(0)
        else:
            self.t1_num.SetValue(1.3)
            self.bat_num.SetValue(1.3)
        self.calibration.update()
        self.update()

    def labelling_changed(self, pasl):
        if pasl:
            self.bat_num.SetValue(0.7)
            self.ie_num.SetValue(0.98)
        else:
            self.bat_num.SetValue(1.3)
            self.ie_num.SetValue(0.85)
