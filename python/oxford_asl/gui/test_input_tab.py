import os
import tempfile

import wx
import nibabel as nib
import numpy as np

from .input_tab import AslInputOptions

class TestInputTab:
    def setup_method(self, method):
        self.app = wx.App(redirect=False)
        self.win = wx.Frame(None)
        #self.app.SetTopWindow(self.win)
        #self.panel = wx.Panel(self.win)
        #self.sizer = wx.BoxSizer(wx.VERTICAL)
        #self.panel.SetSizer(self.sizer)
        self.tab = AslInputOptions(self.win, 0, 1)
        self.tempdir = tempfile.mkdtemp()
        self.asl_data = np.zeros([2, 2, 2, 48], dtype=np.float32)
        self.asl_nii = nib.Nifti1Image(self.asl_data, None)
        self.asl_fname = os.path.join(self.tempdir, "asl.nii.gz")
        self.asl_nii.to_filename(self.asl_fname)
        self.eventloop = wx.EventLoop()

    def process_events(self):
        ea = wx.EventLoopActivator(self.eventloop)
        while self.eventloop.Pending():
            self.eventloop.Dispatch()
        self.eventloop.ProcessIdle()
        del ea

    def test_spld(self):
        self.tab.data_picker.SetPath(self.asl_fname)
        self.tab.update()
        assert(self.tab.nvols == 48)
