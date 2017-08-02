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

def set_localdir(localdir):
    global LOCAL_DIR
    LOCAL_DIR = localdir

class Prog:
    def __init__(self, cmd, localdir=""):
        self.cmd = cmd

    def _find(self):
        """ 
        Find the program, either in the 'local' directory, or in $FSLDEVDIR/bin or $FSLDIR/bin 
        This is called each time the program is run so the caller can control where programs
        are searched for at any time
        """
        if "FSLDIR" in os.environ: 
            fsldir = os.environ["FSLDIR"]
        else:
            fsldir = LOCAL_DIR
        if "FSLDEVDIR" in os.environ: 
            fsldevdir = os.environ["FSLDEVDIR"]
        else:
            fsldevdir = LOCAL_DIR

        local_path = os.path.join(LOCAL_DIR, self.cmd)
        if os.path.isfile(local_path) and os.access(local_path, os.X_OK):
            return local_path
        elif os.path.exists(os.path.join(fsldevdir, "bin/%s" % self.cmd)):
            return os.path.join(fsldevdir, "bin/%s" % self.cmd)
        elif os.path.exists(os.path.join(fsldir, "bin/%s" % self.cmd)):
            return os.path.join(fsldir, "bin/%s" % self.cmd)
        else:
            return self.cmd
    
    def run(self, args):
        """ Run, writing output to stdout and returning retcode """
        cmd = self._find()
        cmd_args = shlex.split(cmd + " " + args)
        out = ""
       #print(" ".join(cmd_args))
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
        self.base = fname.split(".", 1)[0]

        # Try to figure out the extension
        exts = ["", ".nii", ".nii.gz"]
        matches = []
        for ext in exts:
            if os.path.exists("%s%s" % (self.base, ext)):
                matches.append(ext)

        if len(matches) == 0:
            raise RuntimeError("%s file %s not found" % (file_type, fname))
        elif len(matches) > 1:
            raise RuntimeError("%s file %s is ambiguous" % (file_type, fname))
    
        self.ext = matches[0]
        self.full = os.path.abspath("%s%s" % (fname, ext))
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
    Prog("imcp").run("%s %s" % (src, dest))
    
def mkdir(dir, fail_if_exists=False, warn_if_exists=True):
    try:
        os.makedirs(dir)
    except OSError as e:
        if e.errno == errno.EEXIST:
            if fail_if_exists: raise
            elif warn_if_exists: print("WARNING: mkdir - Directory %s already exists" % dir)
