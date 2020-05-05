"""
BASIL GUI for Oxford ASL - Base classes and functions

Copyright (C) 2020 University of Oxford
"""
import os
import functools

import nibabel as nib

class OptionError(RuntimeError):
    """
    Exception thrown because of an invalid option value (i.e.
    an error that the user must resolve)
    """

class OptionComponent():
    """
    Base class for anything which defines / responds to options
    """
    def __init__(self, app):
        self.app = app

    def state_changed(self, _event=None):
        """
        Called when a the state of the component is changed.

        Notify the controlling app which will collect options from all
        component and then call option_changed for each option that has changed.
        """
        self.app.update_options()

    def options(self):
        """
        @return Dictionary of options defined
        """
        return {}

    def option_changed(self, options, key, value):
        """
        Called when an option value is changed

        :param options: Dictionary of options
        :param key: Key of option which has changed
        :param value: New value of option
        """

    def wp_status(self, options, issues):
        """
        :param options: Options dictionary
        :param issues: Sequence of strings to append any deviations from white paper mode in this component
        """

    def check_options(self, options):
        """
        Check the options for consistency and validity

        :param options: Current options dictionary

        :raises OptionError: if there are any problems with the options as defined
        """

    def _check_exists(self, label, fname, can_be_none=True):
        """
        Check if a file exists

        :param label: Human readable label for file (e.g. 'Calibration data')
        :param fname: Filename
        :param can_be_none: If None is an acceptable value

        :raises OptionError: if file does not exist
        """
        if fname is None and can_be_none:
            return
        elif not fname:
            raise OptionError("%s must be specified" % label)
        elif not os.path.exists(fname):
            raise OptionError("%s - no such file or directory" % label)

    @functools.lru_cache(maxsize=128)
    def _check_image(self, label, fname, can_be_none=True):
        """
        :param label: Human readable label for file (e.g. 'Calibration data')
        :param fname: Filename
        :param can_be_none: If None is an acceptable value

        :raises OptionError: if file does not contain an image
        """
        if fname is None and can_be_none:
            return

        self._check_exists(label, fname, can_be_none)
        try:
            nii = nib.load(fname)
        except nib.filebasedimages.ImageFileError:
            raise OptionError("%s - failed to load file - is this a valid image file?" % label)

        try:
            return nii.get_data().shape
        except:
            raise OptionError("%s - failed to read data shape - check image file is not corrupted")

# The options we need to pass to oxford_asl for various data orderings
ORDER_IBF_IAF = {
    "trp"    : ("rpt", "diff"),
    "trp,tc" : ("rpt", "tcb"),
    "trp,ct" : ("rpt", "ctb"),
    "rtp"    : ("tis", "diff"),
    "rtp,tc" : ("tis", "tcb"),
    "rtp,ct" : ("tis", "ctb"),
    "ptr,tc" : ("rpt", "tc"),
    "ptr,ct" : ("rpt", "ct"),
    "prt,tc" : ("tis", "tc"),
    "prt,ct" : ("tis", "ct")
}

def get_ibf_iaf(order, ntc, tagfirst):
    """
    Get the oxford_asl IAF/IBF options corresponding to the
    specified data type and ordering
    """
    if ntc == 2:
        if tagfirst:
            order += ",tc"
        else: order += ",ct"
    if order not in ORDER_IBF_IAF:
        raise OptionError("This data ordering is not supported by ASL_FILE")
    else:
        return ORDER_IBF_IAF[order]

def get_order_ntc_tagfirst(ibf, iaf):
    """
    Get the data ordering, number of TC images and tag-first flag
    for a given set of oxford_asl IBF/IAF options
    """
    order = None
    for order_str, ibf_iaf in ORDER_IBF_IAF.items():
        if ibf == ibf_iaf[0] and iaf == ibf_iaf[1]:
            parts = order_str.split(",")
            order = parts[0]
            if len(parts) == 1:
                ntc = 1
                tagfirst = False
            else:
                ntc = 2
                tagfirst = parts[1] == "tc"
            break

    if order is None:
        raise RuntimeError("Invalid IAF/IBF: %s, %s" % (iaf, ibf))

    return order, ntc, tagfirst

def get_nvols(fname):
    """
    Get the number of volumes in a Nifti data set

    Does not throw an exception on error - use _check_image for that

    :param fname: File name
    :return Number of volumes or -1 if could not read the file for any reason
    """
    try:
        shape = nib.load(fname).get_data().shape
        if len(shape) == 4:
            return shape[3]
        else:
            return 1
    except (FileNotFoundError, nib.filebasedimages.ImageFileError):
        return -1
