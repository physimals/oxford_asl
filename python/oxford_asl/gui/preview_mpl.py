"""
BASIL GUI for Oxford ASL - Preview widgets

Copyright (C) 2020 University of Oxford
"""
import os
import shutil
import tempfile

import numpy as np
import nibabel as nib

import wx

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure

from . import OptionComponent, get_order_ntc_tagfirst
from .cmdline import FslCmd
from .preview_structure import DataStructurePreview

class PreviewPanel(wx.Panel, OptionComponent):
    """
    Panel providing a simple perfusion weighted image preview for the output of ASL_FILE.

    Used so user can check their choice of data grouping/ordering looks right
    """
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, size=wx.Size(300, 600))
        self._options = {"i" : None}
        self._preview_data = None
        self._slice = -1
        self._nslices = 1
        self._view_slice = 0
        self._fig = Figure(figsize=(3.5, 3.5), dpi=100, facecolor='black')
        self._axes = self._fig.add_subplot(111, facecolor='black')
        self._axes.get_xaxis().set_ticklabels([])
        self._axes.get_yaxis().set_ticklabels([])
        self._canvas = FigureCanvas(self, -1, self._fig)
        self._canvas.mpl_connect('scroll_event', self._scroll)
        self._canvas.mpl_connect('button_press_event', self._view_change)
        self._sizer = wx.BoxSizer(wx.VERTICAL)
        font = self.GetFont()
        font.SetWeight(wx.BOLD)
        text = wx.StaticText(self, label="Data preview - perfusion weighted image")
        text.SetFont(font)
        self._sizer.AddSpacer(10)
        self._sizer.Add(text, 0)
        self._sizer.Add(self._canvas, 2, border=5, flag=wx.EXPAND | wx.ALL)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(wx.StaticText(self, label="Use scroll wheel to change slice, double click to change view"), 0, flag=wx.ALIGN_CENTRE_VERTICAL)
        self.update_btn = wx.Button(self, label="Update")
        self.update_btn.Bind(wx.EVT_BUTTON, self._update_clicked)
        hbox.Add(self.update_btn)
        self._sizer.Add(hbox)

        self._sizer.AddSpacer(10)
        text = wx.StaticText(self, label="Data order preview")
        text.SetFont(font)
        self._sizer.Add(text, 0)
        self._data_struc_preview = DataStructurePreview(self, 1, [1], True, "trp", True)
        self._sizer.Add(self._data_struc_preview, 2, wx.EXPAND)
        self.SetSizer(self._sizer)

    def option_changed(self, options, key, value):
        self._options = options

        order, _ntc, _tagfirst = get_order_ntc_tagfirst(options["ibf"], options["iaf"])
        self._data_struc_preview.ntis = options["_ntis"]
        self._data_struc_preview.repeats = options["rpts"]
        self._data_struc_preview.ntc = options["_ntc"]
        self._data_struc_preview.tagfirst = options["iaf"] == "tc"
        self._data_struc_preview.order = order
        self._data_struc_preview.tis_name = "PLDs" if options["casl"] else "TIs"
        self._data_struc_preview.Refresh()

    def _update_clicked(self, _evt):
        """
        Update the preview. This is called explicitly when the user clicks the update
        button as it involves calling ASL_FILE and may be slow
        """
        self._preview_data = None
        infile = self._options["i"]
        if not infile:
            # Don't bother if we have not input file yet!
            return

        tempdir = tempfile.mkdtemp()
        try:
            meanfile = "%s/mean" % tempdir
            cmd = FslCmd("asl_file")
            cmd.add_arg('--data="%s"' % infile)
            cmd.add_arg("--ntis=%i" % self._options["_ntis"])
            cmd.add_arg('--mean="%s"' % meanfile)
            cmd.add_arg("--iaf=%s" % self._options["iaf"])
            cmd.add_arg("--ibf=%s" % self._options["ibf"])
            cmd.run()

            for ext in (".nii", ".nii.gz"):
                if os.path.exists(meanfile + ext):
                    self._preview_data = nib.load(meanfile + ext).get_data()

            if self._preview_data is None:
                raise RuntimeError("Could not load output from asl_file")

            # If multi-TI data, take mean over volumes
            if self._preview_data.ndim == 4:
                self._preview_data = np.mean(self._preview_data, axis=3)

            self._view_slice = 0
            self._init_view()
            self._redraw()
        finally:
            shutil.rmtree(tempdir)

    def _init_view(self):
        self._nslices = self._preview_data.shape[2-self._view_slice]
        self._slice = int(self._nslices / 2)
        self._redraw()

    def _redraw(self):
        """
        Redraw the preview image
        """
        self._axes.clear()
        if self._preview_data is None:
            return

        if self._view_slice == 0:
            slicedata = self._preview_data[:, :, self._slice]
        elif self._view_slice == 1:
            slicedata = self._preview_data[:, self._slice, :]
        else:
            slicedata = self._preview_data[self._slice, :, :]

        img = self._axes.imshow(slicedata.T, interpolation="nearest", vmin=slicedata.min(), vmax=slicedata.max())
        self._axes.set_ylim(self._axes.get_ylim()[::-1])
        img.set_cmap("gray")

        # This rubbish is required to make the display visible without resizing
        self.Layout()
        self.GetParent().Layout()

    def _view_change(self, event):
        """
        Called on mouse click event. Double click changes the view direction and redraws
        """
        if self._preview_data is None:
            return
        if event.dblclick:
            self._view_slice = (self._view_slice + 1) % 3
            self._init_view()
            self._redraw()

    def _scroll(self, event):
        """
        Called on mouse scroll wheel to move through the slices in the current view
        """
        if event.button == "up":
            if self._slice != self._nslices-1:
                self._slice += 1
        else:
            if self._slice != 0:
                self._slice -= 1
        self._redraw()
