#!/usr/bin/env python
"""
Simple wxpython based GUI front-end to OXFORD_ASL command line tool

Currently this does not use any of the FSL python libraries. Possible improvements would be:
  - Use props library to hold run options and use the signalling mechanisms to communicate
    values. The built-in widget builder does not seem to be flexible enough however.
  - Use fsleyes embedded widget as the preview for a nicer, more interactive data preview

Requirements:
  - wxpython 
  - matplotlib
  - numpy
  - nibabel
"""

import sys
import os
import subprocess
import shlex

import wx
import wx.grid

from .analysis_tab import AslAnalysis
from .structure_tab import StructureTab
from .calib_tab import AslCalibration
from .input_tab import AslInputOptions
from .dist_corr_tab import AslDistCorr

from .widgets import PreviewPanel
from .run_box import AslRun

class AslGui(wx.Frame):
    """
    Main GUI window
    """

    def __init__(self):
        wx.Frame.__init__(self, None, title="Basil", size=(1200, 700), style=wx.DEFAULT_FRAME_STYLE)
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
        hsizer.Add(notebook, 1, wx.ALL|wx.EXPAND, 5)
        self.preview = PreviewPanel(hpanel)
        hsizer.Add(self.preview, 1, wx.EXPAND)
        hpanel.SetSizer(hsizer)
        main_vsizer.Add(hpanel, 2, wx.EXPAND)

        self.run_panel = wx.Panel(main_panel)
        runsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.run_label = wx.StaticText(self.run_panel, label="Unchecked")
        self.run_label.SetFont(wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        runsizer.Add(self.run_label, 1, wx.EXPAND)
        self.run_btn = wx.Button(self.run_panel, label="Run")
        runsizer.Add(self.run_btn, 0, wx.ALIGN_CENTER_VERTICAL)
        self.run_panel.SetSizer(runsizer)
        main_vsizer.Add(self.run_panel, 0, wx.EXPAND)
        self.run_panel.Layout()
        
        main_panel.SetSizer(main_vsizer)
        
        self.run = AslRun(self, self.run_btn, self.run_label)
        setattr(self.run, "preview", self.preview)
        tab_cls = [AslInputOptions, StructureTab, AslCalibration, AslDistCorr, AslAnalysis,]
        tabs = [cls(notebook, idx, len(tab_cls)) for idx, cls in enumerate(tab_cls)]
        
        for idx, tab in enumerate(tabs):
            notebook.AddPage(tab, tab.title)
            setattr(tab, "run", self.run)
            setattr(tab, "preview", self.preview)
            setattr(self.run, tab.name, tab)
            setattr(self.preview, tab.name, tab)
            for tab2 in tabs:
                if tab != tab2: setattr(tab, tab2.name, tab2)
            tab.update()

        self.Layout()

def main():
    app = wx.App(redirect=False)
    top = AslGui()
    top.Show()
    app.MainLoop()

if __name__ == '__main__':
    main()
