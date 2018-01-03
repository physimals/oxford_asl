#!/bin/env python
#
# Build script for CMAKE module
#
# Assumes FSLDIR is set. If FSLDEVDIR is also set, this will be used as the installation
# prefix
#
# Note that arch is only currently used for Windows to set the correct VC compiler architecture

import os, sys
import shutil
import traceback
import stat

def remove_readonly(func, path, excinfo):
    os.chmod(path, stat.S_IWRITE)
    func(path)

def rmdir(d):
    try:
        shutil.rmtree(d, onerror=remove_readonly)
    except:
        print("Error removing %s" % d)
        traceback.print_exc(limit=0)

win = sys.platform.startswith("win")
osx = sys.platform.startswith("darwin")

if len(sys.argv) < 3:
    print("Usage: build.py <type> <arch> [--install]")
    sys.exit(1)

build_type = sys.argv[1]
arch = sys.argv[2]
install = "--install" in sys.argv

rootdir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
builddir = os.path.join(rootdir, "build_%s" % build_type)
rmdir(builddir)
os.makedirs(builddir)
cwd = os.getcwd()
os.chdir(builddir)

cmake_opts = "-DCMAKE_BUILD_TYPE=%s" % build_type
if install:
    installdir = os.environ.get("FSLDEVDIR", os.environ["FSLDIR"])  
    cmake_opts += ' -DCMAKE_INSTALL_PREFIX="%s"' % installdir

if win:
    if "VCINSTALLDIR" not in os.environ:
        print("You must run this script from the Visual Studio tools command line")
        sys.exit(1)
    os.system('"%s/vcvarsall" %s' % (os.environ["VCINSTALLDIR"], arch))
    cmake_opts += ' -G "NMake Makefiles"'
    make = "nmake"
else:
    make = "make"

os.system("cmake .. %s" % cmake_opts)
os.system(make)
if install:
    os.system("%s install" % make)
else:
    os.system(make)
    
os.chdir(cwd)

