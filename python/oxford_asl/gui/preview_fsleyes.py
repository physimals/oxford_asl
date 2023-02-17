"""
BASIL GUI for Oxford ASL - Preview widgets

Copyright (C) 2020 University of Oxford
"""
import os
import shutil
import tempfile

import numpy as np
import nibabel as nib

import fsleyes
import fsleyes.overlay as fsloverlay
import fsleyes.displaycontext as fsldc
import fsleyes.views.orthopanel as orthopanel
import fsleyes.profiles as profiles
import fsleyes.colourmaps as colourmaps

from fsl.utils.platform import platform as fslplatform
from fsl.utils import idle
import fsl.data.image as fslimage

import wx

from . import OptionComponent, get_order_ntc_tagfirst
from .cmdline import FslCmd
from .preview_structure import DataStructurePreview

def fsleyes_embed(parent=None, make_fsleyesframe=True, **kwargs):
    """Initialise FSLeyes and create a :class:`.FSLeyesFrame`, when
    running within another application.
    .. note:: If a ``wx.App`` does not exist, one is created.
    :arg parent: ``wx`` parent object
    :arg make_fsleyesframe: bool, default is True to make a new :class:`.FSLeyesFrame`
    :returns:    A tuple containing:
                    - The :class:`.OverlayList`
                    - The master :class:`.DisplayContext`
                    - The :class:`.FSLeyesFrame` or None if make_fsleyesframe=False
    All other arguments are passed to :meth:`.FSLeyesFrame.__init__`.
    """

    import fsleyes_props          as props
    import fsleyes.gl             as fslgl
    import fsleyes.frame          as fslframe
    import fsleyes.overlay        as fsloverlay
    import fsleyes.displaycontext as fsldc

    app    = wx.GetApp()
    ownapp = app is None
    if ownapp:
        app = FSLeyesApp()

    fsleyes.initialise()
    colourmaps.init()
    props.initGUI()

    called = [False]
    ret    = [None]

    def until():
        return called[0]

    def ready():
        frame = None
        fslgl.bootstrap()

        overlayList = fsloverlay.OverlayList()
        displayCtx  = fsldc.DisplayContext(overlayList)
        if make_fsleyesframe:
            frame       = fslframe.FSLeyesFrame(
                parent, overlayList, displayCtx, **kwargs)

        if ownapp:
            app.SetOverlayListAndDisplayContext(overlayList, displayCtx)
            # Keep a ref to prevent the app from being GC'd
            if make_fsleyesframe:
                frame._embed_app = app

        called[0] = True
        ret[0]    = (overlayList, displayCtx, frame)

    fslgl.getGLContext(ready=ready, raiseErrors=True)
    idle.block(10, until=until)

    if ret[0] is None:
        raise RuntimeError('Failed to start FSLeyes')
    return ret[0]

class PreviewPanel(wx.Panel, OptionComponent):
    """
    Panel providing a simple perfusion weighted image preview for the output of ASL_FILE.

    Used so user can check their choice of data grouping/ordering looks right
    """
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, size=wx.Size(300, 600))
        self._options = {"i" : None}
        self._preview_data = None
        self._sizer = wx.BoxSizer(wx.VERTICAL)
        font = self.GetFont()
        font.SetWeight(wx.BOLD)
        text = wx.StaticText(self, label="Data preview - perfusion weighted image")
        text.SetFont(font)
        self._sizer.AddSpacer(10)
        self._sizer.Add(text, 0)
        self.overlayList, masterDisplayCtx, self._fsleyes_frame = fsleyes_embed(None, make_fsleyesframe=False)
        self.displayCtx = fsldc.DisplayContext(self.overlayList, parent=masterDisplayCtx)
        self.op = orthopanel.OrthoPanel(
                self,
                self.overlayList,
                self.displayCtx,
                None)
        self.op.SetMinSize((-1, 300))
        self.op.Show()
        self._sizer.Add(self.op, proportion=5, flag=wx.EXPAND | wx.ALL, border=5)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.update_btn = wx.Button(self, label="Update")
        self.update_btn.Bind(wx.EVT_BUTTON, self._update_clicked)
        hbox.Add(self.update_btn)
        self._sizer.Add(hbox)

        self._sizer.AddSpacer(10)
        text = wx.StaticText(self, label="Data order preview")
        text.SetFont(font)
        self._sizer.Add(text, 0)
        self._data_struc_preview = DataStructurePreview(self, 1, [1], True, "trp", True)
        self._sizer.Add(self._data_struc_preview, 1, wx.EXPAND)
        self.SetSizer(self._sizer)

    def OnClose(self, evt=None):
        self._fsleyes_frame._embed_app.Exit()

    def option_changed(self, options, key, value):
        self._options = options

        order, _ntc, _tagfirst = get_order_ntc_tagfirst(options["ibf"], options["iaf"])
        self._data_struc_preview.ntis = options["_ntis"]
        self._data_struc_preview.repeats = options["rpts"]
        self._data_struc_preview.ntc = options["_ntc"]
        self._data_struc_preview.tagfirst = options["iaf"] == "tc"
        self._data_struc_preview.order = order
        self._data_struc_preview.tis_name = "PLD" if options["casl"] else "TI"
        self._data_struc_preview.Refresh()

    def _update_clicked(self, _evt):
        """
        Update the preview. This is called explicitly when the user clicks the update
        button as it involves calling ASL_FILE and may be slow
        """
        self._preview_data = None
        self.overlayList.clear()
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
            if self._options["iaf"] != "diff":
                cmd.add_arg("--diff")
            cmd.run()

            for ext in (".nii", ".nii.gz"):
                if os.path.exists(meanfile + ext):
                    self._preview_data = nib.load(meanfile + ext).get_fdata()

            if self._preview_data is None:
                raise RuntimeError("Could not load output from asl_file")

            # If multi-TI data, take mean over volumes
            if self._preview_data.ndim == 4:
                self._preview_data = np.mean(self._preview_data, axis=3)

            self._view_slice = 0
            self._redraw()
        finally:
            shutil.rmtree(tempdir)

    def _redraw(self):
        """
        Redraw the preview image
        """
        if self._preview_data is None:
            return

        img = fslimage.Image(self._preview_data)
        self.overlayList.append(img)
