#!/usr/bin/env python
#
# FSL helper functions - named fslhelpers so as not to clash with 'official' FSL python
# modules

import os, sys
import glob
import shutil
import shlex
import subprocess
import math
import errno

import nibabel as nib
import numpy as np

# Where to look for 'local' programs - default to dir of calling program
# This enables users to override standard FSL programs by putting their
# own copies in the same directory as the script they are calling
LOCAL_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))
#print("FSLHELPERS: using local binaries dir: %s" % LOCAL_DIR)

ECHO = False

def set_localdir(localdir):
    global LOCAL_DIR
    LOCAL_DIR = localdir

def set_echo(echo):
    global ECHO
    ECHO = echo

class Prog:
    def __init__(self, cmd):
        self.cmd = cmd

    def _find(self):
        """ 
        Find the program, either in the 'local' directory, or in $FSLDEVDIR/bin or $FSLDIR/bin 
        This is called each time the program is run so the caller can control where programs
        are searched for at any time
        """
        dirs = [
            LOCAL_DIR,
            os.path.join(os.environ.get("FSLDEVDIR", LOCAL_DIR), "bin"),
            os.path.join(os.environ.get("FSLDIR", LOCAL_DIR), "bin"), 
        ]

        for d in dirs:
            ex = os.path.join(d, self.cmd)
            if os.path.isfile(ex) and os.access(ex, os.X_OK):
                return ex
        
        return self.cmd
    
    def run(self, args):
        return self(args)

    def __call__(self, args):
        """ Run, writing output to stdout and returning retcode """
        cmd = self._find()
        cmd_args = shlex.split(cmd + " " + args)
        out = ""
        global ECHO
        if ECHO:
            print(" ".join(cmd_args))
        p = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        while 1:
            retcode = p.poll() #returns None while subprocess is running
            line = p.stdout.readline()
            out += line
            if retcode is not None: break
        if retcode != 0:
            raise RuntimeError(out)
        else:
            return out

class Image: 
    def __init__(self, fname, file_type="Image"):
        if fname is None or fname == "":
            raise RuntimeError("%s file not specified" % file_type)

        self.file_type = file_type
        d, name = os.path.split(fname)
        self.base = name.split(".", 1)[0]
        self.noext = os.path.join(os.path.abspath(d), self.base)

        # Try to figure out the extension
        exts = ["", ".nii", ".nii.gz"]
        matches = []
        for ext in exts:
            if os.path.exists("%s%s" % (self.noext, ext)):
                matches.append(ext)

        if len(matches) == 0:
            raise RuntimeError("%s file %s not found" % (file_type, fname))
        elif len(matches) > 1:
            raise RuntimeError("%s file %s is ambiguous" % (file_type, fname))
    
        self.ext = matches[0]
        self.full = os.path.abspath("%s%s" % (self.noext, ext))
        self.nii = nib.load(self.full)
        self.shape = self.nii.shape
         
    def data(self):
        return self.nii.get_data()

    def new_nifti(self, data=None):
        """ Return new Nifti oject, taking header info from this image """
        if data is None: data=self.data()
        nii = nib.Nifti1Image(data, self.nii.header.get_best_affine())
        nii.update_header()
        return nii

    def check_shape(self, shape=None):
        if len(self.shape) != len(shape):
            raise RuntimeError("%s: expected %i dims, got %i" % (self.file_type, len(shape), len(self.shape)))
        if self.shape != shape:
            raise RuntimeError("%s: shape (%s) does not match (%s)" % (self.file_type, str(self.shape), str(shape)))

maths = Prog("fslmaths")
roi = Prog("fslroi")
stats = Prog("fslstats")
merge = Prog("fslmerge")
bet = Prog("bet")
flirt = Prog("flirt")
fast = Prog("fast")

def imcp(src, dest):
    Prog("imcp")("%s %s" % (src, dest))
    
def mkdir(dir, fail_if_exists=False, warn_if_exists=True):
    try:
        os.makedirs(dir)
    except OSError as e:
        if e.errno == errno.EEXIST:
            if fail_if_exists: raise
            elif warn_if_exists: print("WARNING: mkdir - Directory %s already exists" % dir)

def tmpdir(suffix, debug=False):
    if debug:
        tmpdir = os.path.join(os.getcwd(), "tmp_%s" % suffix)
        fsl.mkdir(tmpdir)
    else:
        tmpdir = tempfile.mkdtemp("_%s" % suffix)
    print("Using temporary dir: %s" % tmpdir)
    return tmpdir