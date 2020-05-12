"""
BASIL GUI for Oxford ASL - Code to run oxford_asl script and display output

Copyright (C) 2020 University of Oxford
"""

import os

import wx
from wx.lib.pubsub import pub

from . import OptionComponent
from .cmdline import CmdRunner, Mkdir, FslCmd

class OxfordAslRunner(wx.Frame, OptionComponent):
    """
    Runs oxford_asl (and any other required commands) based on configured options
    """

    def __init__(self, parent):
        wx.Frame.__init__(self, parent, title="Run", size=(600, 400), style=wx.DEFAULT_FRAME_STYLE)
        self._run_sequence = []

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.output_text = wx.TextCtrl(self, style=wx.TE_READONLY | wx.TE_MULTILINE)
        font = wx.Font(8, wx.TELETYPE, wx.NORMAL, wx.NORMAL)
        self.output_text.SetFont(font)
        self.sizer.Add(self.output_text, 1, flag=wx.EXPAND)

        self.SetSizer(self.sizer)
        self.Bind(wx.EVT_CLOSE, self._close)
        pub.subscribe(self._write_output, "run_stdout")
        pub.subscribe(self._finished, "run_finished")

    def _write_output(self, line):
        self.output_text.AppendText(line)

    def _close(self, _):
        self.Hide()

    def _finished(self, retcode):
        if retcode != 0:
            self._write_output("\nWARNING: command failed\n")

    def run(self):
        self.Show()
        self.Raise()
        self.output_text.Clear()
        runner = CmdRunner(self._run_sequence)
        runner.start()

    def option_changed(self, options, key, value):
        """
        Get the sequence of commands for the selected options, throwing exception
        if any problems are found (e.g. files don't exist, mandatory options not specified)

        Exception text is reported by the GUI
        """
        self._run_sequence = []
        options = dict(options) # We may change it

        # Create output directory if required
        outdir = options["o"]
        self._run_sequence.append(Mkdir(outdir))

        # Run FSL_ANAT if required
        fslanat = options.get("fslanat", None)
        if fslanat is not None and not os.path.exists(fslanat):
            fslanat_options = {
                "i" : options.pop("s"),
                "o" : "%s/struc" % outdir,
            }
            fslanat_cmd = FslCmd("fsl_anat", fslanat_options)
            self._run_sequence.append(fslanat_cmd)

        # Finally run OXFORD_ASL command
        cmd = FslCmd("oxford_asl", options)
        self._run_sequence.append(cmd)
