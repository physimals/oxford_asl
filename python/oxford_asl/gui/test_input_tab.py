import os
import tempfile

import wx
import nibabel as nib
import numpy as np

from .input_tab import AslInputOptions
from .main import AslGui

class MockApp(wx.App):

    def update_options(self):
        pass

class TestInputTab:
    def setup_method(self, method):
        app = wx.App(redirect=False)
        self.app = AslGui()
        self.win = wx.Frame(None)
        #self.app.SetTopWindow(self.win)
        #self.panel = wx.Panel(self.win)
        #self.sizer = wx.BoxSizer(wx.VERTICAL)
        #self.panel.SetSizer(self.sizer)
        self.tab = AslInputOptions(self.app, self.app, 0, 1)
        self.tempdir = tempfile.mkdtemp()
        self.asl_data = np.zeros([2, 2, 2, 48], dtype=np.float32)
        self.asl_nii = nib.Nifti1Image(self.asl_data, None)
        self.asl_fname = os.path.join(self.tempdir, "asl.nii.gz")
        self.asl_nii.to_filename(self.asl_fname)
        #self.eventloop = wx.EventLoop()

    def test_spld(self):
        self.tab.data_picker.SetPath(self.asl_fname)
        self.tab.state_changed()
        options = self.tab.options()
        assert(options["_nvols"] == 48)
