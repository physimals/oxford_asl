#!/usr/bin/env python
"""
BASIL GUI for Oxford ASL - Main window

Currently this does not use any of the FSL python libraries. Possible improvements would be:
  - Use props library to hold run options and use the signalling mechanisms to communicate
    values. The built-in widget builder does not seem to be flexible enough however.
  - Use fsleyes embedded widget as the preview for a nicer, more interactive data preview

Requirements:
  - wxpython
  - fsleyes (or matplotlib using --matplotlib option)
  - numpy
  - nibabel
"""
import sys
import os
import traceback

import wx
import wx.grid

from . import OptionError

from .widgets import WhitePaperCompatibility
from .analysis_tab import AnalysisTab
from .structure_tab import StructureTab
from .calib_tab import CalibTab
from .input_tab import AslInputOptions
from .dist_corr_tab import DistCorrTab

from .runner import OxfordAslRunner

class AslGui(wx.Frame):
    """
    Main GUI window
    """

    def __init__(self):
        wx.Frame.__init__(self, None, title="Basil", size=(1200, 750), style=wx.DEFAULT_FRAME_STYLE)
        icon_fname = os.path.join(os.path.abspath(os.path.dirname(__file__)), "basil.png")
        self.SetIcon(wx.Icon(icon_fname))
        self._options = {}
        main_panel = wx.Panel(self)
        main_vsizer = wx.BoxSizer(wx.VERTICAL)

        banner = wx.Panel(main_panel, size=(-1, 80))
        banner.SetBackgroundColour((54, 122, 157))
        banner_fname = os.path.join(os.path.abspath(os.path.dirname(__file__)), "banner.png")
        wx.StaticBitmap(banner, -1, wx.Bitmap(banner_fname, wx.BITMAP_TYPE_ANY))
        main_vsizer.Add(banner, 0, wx.EXPAND)

        hpanel = wx.Panel(main_panel)
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        notebook = wx.Notebook(hpanel, id=wx.ID_ANY, style=wx.BK_DEFAULT)
        hsizer.Add(notebook, 0, wx.ALL, 5)

        if "--matplotlib" in sys.argv:
            from .preview_mpl import PreviewPanel
        else:
            from .preview_fsleyes import PreviewPanel
        preview = PreviewPanel(hpanel)
        hsizer.Add(preview, 1, wx.EXPAND)
        hpanel.SetSizer(hsizer)
        main_vsizer.Add(hpanel, 1, wx.EXPAND)

        line = wx.StaticLine(main_panel, style=wx.LI_HORIZONTAL)
        main_vsizer.Add(line, 0, wx.EXPAND)

        bottom_panel = wx.Panel(main_panel)
        bottom_sizer = wx.BoxSizer(wx.HORIZONTAL)
        bottom_panel.SetSizer(bottom_sizer)
        
        self.run_btn = wx.Button(bottom_panel, label="Run")
        bottom_sizer.Add(self.run_btn, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        self.run_label = wx.StaticText(bottom_panel, label="Unchecked")
        self.run_label.SetFont(wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        bottom_sizer.Add(self.run_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        bottom_sizer.AddStretchSpacer(1)
        self.wpcompat = WhitePaperCompatibility(self, bottom_panel)
        bottom_sizer.Add(self.wpcompat, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        
        main_vsizer.Add(bottom_panel, 0, wx.EXPAND)
        main_panel.SetSizer(main_vsizer)
        runner = OxfordAslRunner(self, self.run_btn, self.run_label)

        self.widgets = [preview, runner, self.wpcompat]
        tab_cls = [AslInputOptions, StructureTab, CalibTab, DistCorrTab, AnalysisTab]
        for idx, cls in enumerate(tab_cls):
            tab = cls(self, notebook, idx, len(tab_cls))
            self.widgets.append(tab)
            notebook.AddPage(tab, tab.title)

        self.update_options()
        self.SetMinSize(self.GetSize())
        #self.SetMaxSize(self.GetSize())
        self.Bind(wx.EVT_CLOSE, self._close)

    def _close(self, evt=None):
        # For some reason with the fsleyes preview the app doesn't exit when then window
        # sis closed. So we will force it to do so...
        self.Destroy()
        sys.exit(0)

    def update_options(self):
        """
        Get the sequence of commands and enable the run button if options are valid. Otherwise
        display the first error in the status label
        """
        try:
            while 1:
                options = {}
                for widget in self.widgets:
                    options.update(widget.options())

                key_diffs = set(options.keys()) ^ set(self._options.keys())
                for key, value in options.items():
                    if value != self._options.get(key, None):
                        key_diffs.add(key)
                if not key_diffs:
                    break

                for diff_key in key_diffs:
                    for widget in self.widgets:
                        widget.option_changed(options, diff_key, options.get(diff_key, None))

                self._options = options

            try:
                for widget in self.widgets:
                    widget.check_options(self._options)

                self.run_label.SetForegroundColour(wx.Colour(0, 128, 0))
                self.run_label.SetLabel("Ready to Go")
                self.run_btn.Enable(True)
            except OptionError as exc:
                self.run_btn.Enable(False)
                self.run_label.SetForegroundColour(wx.Colour(255, 0, 0))
                self.run_label.SetLabel(str(exc))
        except Exception:
            # Any uncaught exception is a program bug - report it to STDERR
            self.run_btn.Enable(False)
            self.run_label.SetForegroundColour(wx.Colour(255, 0, 0))
            self.run_label.SetLabel("Unexpected error - see console and report as a bug")
            traceback.print_exc()

def main():
    """
    GUI entry point
    """
    app = wx.App(redirect=False)
    top = AslGui()
    top.Show()
    app.MainLoop()

if __name__ == '__main__':
    main()
