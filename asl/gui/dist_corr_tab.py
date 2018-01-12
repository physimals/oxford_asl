import wx
import wx.grid

from .widgets import TabPage

class AslDistCorr(TabPage):
    """
    Tab page containing distortion correction options
    """
    FIELDMAP = 0
    CALIB_IMAGE = 1

    def __init__(self, parent, idx, n):
        TabPage.__init__(self, parent, "Distortion Correction", idx, n, name="distcorr")

        self.distcorr_choices = ["Fieldmap", "Calibration image"]

        self.section("Distortion Correction")

        # Calibration image options
        self.distcorr_cb = wx.CheckBox(self, label="Apply distortion correction")
        self.distcorr_cb.Bind(wx.EVT_CHECKBOX, self.update)

        self.distcorr_ch = wx.Choice(self, choices=self.distcorr_choices[:1])
        self.distcorr_ch.SetSelection(self.FIELDMAP)
        self.distcorr_ch.Bind(wx.EVT_CHOICE, self.update)
        self.pack("", self.distcorr_cb, self.distcorr_ch, enable=False)

        # Calib image options
        self.section("Calibration Image Mode")
        self.calib_picker = self.file_picker("Phase-encode-reversed calibration image")
        
        # Fieldmap options
        self.section("Fieldmap Mode")
        self.fmap_picker = self.file_picker("Fieldmap image (in rad/s)")
        self.fmap_mag_picker = self.file_picker("Fieldmap magnitude image")
        self.fmap_be_picker = self.file_picker("Brain-extracted magnitude image", optional=True)
        
        # General options
        self.section("General")
        self.echosp_num = self.number("Effective EPI echo spacing", min=0, max=10)
        self.pedir_ch = self.choice("Phase encoding direction", choices=["x", "y", "z", "-x", "-y", "-z"])
        
        self.sizer.AddGrowableCol(1, 1)
        #sizer.AddGrowableRow(5, 1)
        self.SetSizer(self.sizer)
        self.next_prev()

    def distcorr(self): return self.distcorr_cb.IsChecked()
    def distcorr_type(self): return self.distcorr_ch.GetSelection()

    def calib(self): return self.calib_picker.GetPath()

    def fmap(self): return self.fmap_picker.GetPath()
    def fmap_mag(self): return self.fmap_mag_picker.GetPath()
    def fmap_mag_be(self): return self.fmap_be_picker.GetPath()

    def echosp(self): return self.echosp_num.GetValue()
    def pedir(self): return self.pedir_ch.GetStringSelection()

    def update(self, event=None):
        self.distcorr_ch.Enable(self.distcorr())

        cal = self.distcorr() and self.distcorr_type() == self.CALIB_IMAGE
        self.calib_picker.Enable(cal)

        fmap = self.distcorr() and self.distcorr_type() == self.FIELDMAP
        self.fmap_picker.Enable(fmap)
        self.fmap_mag_picker.Enable(fmap)
        self.fmap_be_picker.Enable(fmap)

        self.pedir_ch.Enable(self.distcorr())
        self.echosp_num.Enable(self.distcorr())

        TabPage.update(self)
        
    def calib_changed(self, enabled):
        """ If calibration enabled, add the calibration image option for distortion correction"""
        sel = self.distcorr_ch.GetSelection()
        if enabled: 
            choices = self.distcorr_choices
            sel = 1
        else: 
            choices = self.distcorr_choices[:1]
            sel = 0
        self.distcorr_ch.Enable(False)
        self.distcorr_ch.Clear()
        self.distcorr_ch.AppendItems(choices)
        self.distcorr_ch.SetSelection(sel)
        self.update()
