"""
BASIL GUI for Oxford ASL - Code to run the oxford_asl pipeline

Copyright (C) 2020 University of Oxford
"""
import os
import sys
import shlex
import subprocess
from threading import Thread

import wx
from wx.lib.pubsub import pub

class Cmd:
    """ Run a shell command """
    def run(self):
        """ Run the command and return exit status """
        return 0

class Mkdir(Cmd):
    """
    Create directories
    """
    def __init__(self, dirname):
        self.dirname = dirname

    def run(self):
        if not os.path.exists(self.dirname):
            os.makedirs(self.dirname)
        return 0

class FslCmd(Cmd):
    """
    Run an FSL command
    """
    def __init__(self, cmd, options=None):
        script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
        fsldevdir = os.path.join(os.environ.get("FSLDEVDIR", ""), "bin")
        fsldir = os.path.join(os.environ.get("FSLDIR", ""), "bin")

        self.cmd = cmd
        for bindir in (script_dir, fsldevdir, fsldir):
            if os.path.isfile(os.path.join(bindir, cmd)):
                self.cmd = os.path.join(bindir, cmd)
                break

        if options is None:
            options = {}

        for key, value in options.items():
            if key.startswith("_"):
                continue

            if len(key) == 1:
                key = "-%s" % key
            else:
                key = "--%s" % key

            if value is None:
                continue
            elif isinstance(value, float):
                value = "%.5g" % value
            elif isinstance(value, list):
                if len(value) > 0:
                    if isinstance(value[0], int):
                        value = ",".join(["%i" % num for num in value])
                    elif isinstance(value[0], float):
                        value = ",".join(["%.5g" % num for num in value])
                    else:
                        value = ",".join([str(num) for num in value])
            elif isinstance(value, bool):
                # Boolean arguments have no value - just include argument if true
                if value:
                    value = None
                else:
                    continue
            else:
                value = str(value)

            self.add_arg(key, value)

    def add_arg(self, key, value=None):
        """
        Add a command line argument
        """
        if value is not None:
            # Filenames might contain spaces so quote them
            if " " in value:
                value = '"%s"' % value
            self.cmd += " %s=%s" % (key, str(value))
        else:
            self.cmd += " %s" % key

    def run(self):
        self._send_output(self.cmd + "\n")
        args = shlex.split(self.cmd)
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        while 1:
            retcode = process.poll() #returns None while subprocess is running
            line = process.stdout.readline().decode('utf-8')
            self._send_output(line)
            if retcode is not None:
                break
        self._send_output("\nReturn code: %i\n\n" % retcode)
        return retcode

    def _send_output(self, line):
        wx.CallAfter(pub.sendMessage, "run_stdout", line=line)

    def __str__(self):
        return self.cmd

class CmdRunner(Thread):
    """
    Runs a sequence of commands in the background
    """
    def __init__(self, cmds):
        """
        :param cmds: Sequence of ``Cmd`` instances
        """
        Thread.__init__(self)
        self.cmds = cmds

    def run(self):
        ret = -1
        try:
            for cmd in self.cmds:
                ret = -1
                ret = cmd.run()
                if ret != 0:
                    break
        finally:
            wx.CallAfter(pub.sendMessage, "run_finished", retcode=ret)
