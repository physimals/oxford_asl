#!/usr/bin/env python
"""
Generates perfusion stats within various ROIs
"""
import os
import sys
import argparse
import csv
import glob

from fsl.data import atlases
from fsl.data.image import Image
import fsl.wrappers as fsl

import nibabel as nib
import numpy as np
import scipy.stats
from scipy.fftpack import fft, ifft

class ArgumentParser(argparse.ArgumentParser):
    """
    ArgumentParser for program options
    """

    def __init__(self, **kwargs):
        argparse.ArgumentParser.__init__(self, prog="oxford_asl_roi_stats", add_help=False, **kwargs)
        self.add_argument("--oxasl-output", required=True,
                          help="OXFORD_ASL or OXASL native space output directory")
        self.add_argument("--add-arrival", action="store_true", default=False,
                          help="Generate output for arrival time as well as perfusion")
        self.add_argument("--fslanat",
                          help="FSL_ANAT output directory - if not specified --struc --gm-pve --wm-pve and --struc2std must be given")
        self.add_argument("--struc", "-s", help="Structural space reference image - ignored if --fslanat is given")
        self.add_argument("--gm-pve", help="GM PVE, assumed to be in structural space unless --native_pves option is specified - ignored if --fslanat is given")
        self.add_argument("--wm-pve", help="WM PVE, assumed to be in structural space unless --native_pves option is specified - ignored if --fslanat is given")
        self.add_argument("--csf-pve", help="CSF PVE, assumed to be in structural space unless --native_pves option is specified - ignored if --fslanat is given")
        self.add_argument("--native-pves", action='store_true', default=False,
                          help="If specified, it is assumed that the GM and WM PVEs provided are in native asl space - ignored if --fslanat is given")
        self.add_argument("--psf", help="Point-spread function for ROI sets. If specified, PSF will be applied to ROI sets (in native space) and 'fuzzy' mean used")
        self.add_argument("--asl2struc", help="File containing ASL->Structural transformation matrix - if not specified will look in <oxasl_output>/asl2struct.mat")
        self.add_argument("--struc2std", help="Structural -> standard space nonlinear warp map - ignored if --fslanat is given. Only required if --std2struc is not given")
        self.add_argument("--std2struc", help="Standard -> structural space nonlinear warp map - ignored if --fslanat is given. If not specified will be derived by inverting --struc2std warp")
        self.add_argument("--output", "-o", required=True,
                          help="Output directory")
        self.add_argument("--min-nvoxels", default=10, type=int,
                          help="Minimum number of relevant voxels required to report statistics")
        self.add_argument("--gm-thresh", default=0.8, type=float,
                          help="Probability threshold for 'pure' grey matter")
        self.add_argument("--wm-thresh", default=0.9, type=float,
                          help="Probability threshold for 'pure' white matter")
        self.add_argument("--min-gm-thresh", default=0.1, type=float,
                          help="Probability threshold for a voxel to be included in GM stats")
        self.add_argument("--min-wm-thresh", default=0.1, type=float,
                          help="Probability threshold for a voxel to be included in WM stats")
        self.add_argument("--roi-native", nargs="*", default=[],
                          help="Additional ROI as binarised mask in ASL space. The name of the ROI will be the stripped filename. May be specified multiple times")
        self.add_argument("--roi-struct", nargs="*", default=[],
                          help="Additional ROI as binarised mask in structural space. The name of the ROI will be the stripped filename. May be specified multiple times")
        self.add_argument("--add-mni-atlas", nargs="*", default=[],
                          help="Additional regions as labelled image in MNI space. If a single region is contained, the name of the ROI will be the stripped filename, otherwise use --add-mni-atlas-names. May be specified multiple times")
        self.add_argument("--add-mni-atlas-labels", nargs="*", default=[],
                          help="Filename containing names of regions in atlas given in --add-mni-atlas")
        self.add_argument("--add-standard-atlases", action="store_true", default=False,
                          help="Add ROIs from Harvard-Oxford cortical/subcortical atlases")
        self.add_argument("--save-mni-rois", action="store_true", default=False,
                          help="Save ROIs in MNI space")
        self.add_argument("--save-struct-rois", action="store_true", default=False,
                          help="Save ROIs in structural space")
        self.add_argument("--save-native-rois", action="store_true", default=False,
                          help="Save ROIs in native (ASL) space")
        self.add_argument("--output-prefix", help="Prefix for output files",
                          default="region_analysis")

def _transform(img, warp, ref, premat=None, postmat=None, interp="trilinear", paddingsize=1, output_is_roi=False, output_roi_thresh=0.5):
    """
    Transform an image

    :param img: fsl.data.Image containing image to transform
    :param warp: Transformation matrix or warp image
    :param ref:  fsl.data.Image containing reference image
    :param premat: Optional Pre-warp affine transformation matrix
    :param premat: Optional Post-warp affine transformation matrix
    :param interp: Interpolation method
    :param paddingsize: Padding size in pixels
    :param output_is_roi: Output should be binarized as an ROI
    :param output_roi_threshold: Thresholding value for binarizing output ROI

    :return:  fsl.data.Image containing transformed image
    """
    kwargs = {
        "warp" : warp, 
        "premat" : premat, 
        "postmat" : postmat,
        "rel" : True,
        "super" : True,
        "superlevel" : "a",
        "interp" : interp,
        "paddingsize" : paddingsize,
    }
    ret = fsl.applywarp(img, ref, out=fsl.LOAD, **kwargs)["out"]

    if output_is_roi:
        # Binarise mask images
        ret = Image((ret.data > output_roi_thresh).astype(np.int), header=ret.header)
    return ret

def _write_nii(img, fname, header=None, vol=None):
    os.makedirs(os.path.dirname(fname), exist_ok=True)
    if isinstance(img, Image):
        header = img.header
        img = img.data

    if img.ndim == 4 and vol is not None:
        img = img[..., vol]
    elif img.ndim == 3 and vol > 0:
        raise RuntimeError("Tried to save multiple volumes of 3D data")

    nii = nib.Nifti1Image(img, None, header=header)
    nii.update_header()
    nii.to_filename(fname)

def _addvar(f):
    """ Add an unused 'variance' parameter to a function which doesn't use it """
    def _f_with_var(val, var):
        return f(val)
    return _f_with_var

def mean_invvarweighted(val, var):
    """ Inverse variance weighted mean (i.e. precision weighted mean) """
    prec = 1 / var
    prec[~np.isfinite(prec)] = 0
    return np.sum(val * prec) / np.sum(prec)

def i2(val, var, mean_precweighted=None):
    """ I^2 Measure of heterogenaity """
    prec = 1 / var
    prec[~np.isfinite(prec)] = 0
    n = len(val)
    if mean_precweighted is not None:
        mu_bar = mean_precweighted
    else:
        mu_bar = mean_invvarweighted(val, var)

    Q = np.sum(prec * (val - mu_bar)**2)
    if Q == 0:
        i2 = 0
    else:
        i2 = (Q-(n-1))/Q

    # Negative values map to 0 (see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC192859/)
    i2 = max(i2, 0)

    # Return I2 as an integer percentage
    return int(100*i2+0.5)

STATS_FNS = {
    "Mean" : _addvar(np.mean),
    "Std" : _addvar(np.std),
    "Median" : _addvar(np.median),
    "IQR" : _addvar(scipy.stats.iqr),
    "Precision-weighted mean" : mean_invvarweighted,
    "I2" : i2,
}

def get_stats_binary(stats, img, var_img, roi, suffix="", ignore_nan=True, ignore_inf=True, ignore_zerovar=True, min_nvoxels=10, mask=None):
    """
    Get a set of statistics for a 3D image within an roi

    :param img: 3D Numpy array
    :param roi: 3D Numpy array with same dimensions as img and boolean data type
    :param ignore_nan: Voxels with care NaN in img are ignored
    :param ignore_inf: Voxels which are infinite in img are ignored
    :param min_nvoxels: If the number of voxels in the ROI is less than this number
                       (after removing Nan and infinte values) no value will be returned

    :return: Mapping from name of statistic to value. This may be NaN or infinite depending
             on the input arguments. If the number of eligible voxels is less than min_nvoxels,
             None is returned (not NaN or zero).
    """
    if list(img.shape) != list(roi.shape):
        raise ValueError("Image must have same dimensions as ROI")
    if list(var_img.shape) != list(roi.shape):
        raise ValueError("Variance image must have same dimensions as ROI")
    if mask is not None and list(mask.shape) != list(roi.shape):
        raise ValueError("Mask must have same dimensions as ROI")

    if mask is None:
        mask = np.ones(roi.shape, dtype=np.int)
    if ignore_nan:
        mask = np.logical_and(mask, ~np.isnan(img))
    if ignore_inf:
        mask = np.logical_and(mask, np.isfinite(img))
    if ignore_zerovar:
        mask = np.logical_and(mask, var_img > 0)

    # Only take voxels where at least one of the ROIs has non-zero percentage
    effective_roi = np.logical_and(roi, mask)

    sample_data = img[effective_roi]
    sample_var = var_img[effective_roi]
    # Variance should not be zero but sometimes is - maybe masking?
    sample_var[sample_var == 0] = 1e-6
    nvoxels = len(sample_data)
    stats["Nvoxels" + suffix] = nvoxels
    for stat, fn in STATS_FNS.items():
        if nvoxels < min_nvoxels:
            stats[stat + suffix] = None
        else:
            stats[stat + suffix] = fn(sample_data, sample_var)

def standardise(roi_set, mode="expand"):
    """
    Ensure that fuzzy ROI set has PVs that sum to 1 in every voxel

    :param mode: 'expand' adds an ROI to the set containing any non-included PV.
                 'normalise' scales PVs so they sum to 1
    """
    if mode not in ("normalise", "expand"):
        raise ValueError(f"Mode {mode} not supported - expected `normalise` or `expand`")

    # Warning if normalising and ROI only has one label
    if np.shape(roi_set)[1] == 1 and mode == "normalise":
        print("Warning: Only one ROI in set with mode='normalise' - all nonzero voxels will have PV=1")

    # Get total PV in each voxel
    roi_sum = np.sum(roi_set, axis=1, keepdims=True)

    if mode == "normalise":
        roi_set = np.where(roi_sum!=0, roi_set/roi_sum, 0)
    elif mode == "expand":
        # TODO: insert a check to see if any diff_img < 0 and raise a warning
        roi_set = np.concatenate((roi_set, 1 - roi_sum), axis=-1)
    return roi_set

def get_stats_fuzzy(stats, img, var_img, roi_set, suffix="", ignore_nan=True, ignore_inf=True, ignore_zerovar=True, mask=None, pv_threshold=0.):
    """
    Get a set of statistics for a set of 'fuzzy' ROIs

    :param img: 3D Numpy array
    :param roi_set: 4D Numpy array with same dimensions as img and each volume
                    containing partial volumes for each ROI in the set

    :return: Mapping from name of statistic to sequence of values, one for each ROI in the set.
             This may be NaN or infinite depending on the input arguments. 
    """
    roi_shape = list(roi_set.shape)[:3]
    if list(img.shape) != roi_shape:
        raise ValueError("Image must have same dimensions as ROI")
    if list(var_img.shape) != roi_shape:
        raise ValueError("Variance image must have same dimensions as ROI")
    if mask is not None and list(mask.shape) != roi_shape:
        raise ValueError("Mask must have same dimensions as ROI")

    if mask is None:
        mask = np.ones(roi_shape, dtype=np.int)
    if ignore_nan:
        mask = np.logical_and(mask, ~np.isnan(img))
    if ignore_inf:
        mask = np.logical_and(mask, np.isfinite(img))
    if ignore_zerovar:
        mask = np.logical_and(mask, var_img > 0)

    # Only take voxels where at least one of the ROIs has non-zero percentage
    mask = np.logical_and(mask, np.sum(roi_set, axis=3) > pv_threshold)

    # Flatten ROI PVs and data into masked 2D array
    roi_array = roi_set[mask]
    g = img[mask]

    # Standardize ROI set so total PV is 1
    roi_array = standardise(roi_array, mode='expand')

    # Ask Jack about this???
    #if var:
    #    roi_array = np.square(roi_array)

    HT = roi_array.T
    print(f" - Fuzzy ROI set: condition number for transfer matrix (unweighted) = {np.linalg.cond(HT):.2f}.")

    # Calculate roi means by linear regression
    means_lstsq, _res, _rank, _s = np.linalg.lstsq(HT@roi_array, HT@g[..., np.newaxis],
                                                   rcond=None) # None uses future default
                                                               # and silences warning

    # Note that we do not report stats for the 'background' ROI added to ensure total PV of 1
    stats["Nvoxels" + suffix] = [np.count_nonzero(roi_array[:, idx] > pv_threshold) for idx in range(roi_set.shape[-1])]
    stats["Mean" + suffix] = np.atleast_1d(np.squeeze(means_lstsq[:-1]))

    # If variance has been supplied add a precision-weighted mean
    if var_img is not None:
        V_inv = scipy.sparse.diags(1/var_img[mask])
        HT = roi_array.T @ V_inv
        print(f" - Fuzzy ROI set: condition number for transfer matrix (prec-weighted) = {np.linalg.cond(HT):.2f}.")

        # Calculate roi means by linear regression
        means_lstsq, _res, _rank, _s = np.linalg.lstsq(HT@roi_array, HT@g[..., np.newaxis],
                                                       rcond=None) # None uses future default
                                                                   # and silences warning
        stats["Precision-weighted mean" + suffix] = np.atleast_1d(np.squeeze(means_lstsq[:-1]))

def get_stats(roi, data_item, options):
    if "fuzzy_native" in roi:
        get_stats_fuzzy(roi["stats"], data_item["f"].data, data_item["var"].data, roi["fuzzy_native"], mask=data_item["mask"])
    elif "mask_native" in roi:
        #print(" - Binarised ROI: %s" % roi["name"])
        get_stats_binary(roi["stats"], data_item["f"].data, data_item["var"].data, roi["mask_native"], mask=data_item["mask"], min_nvoxels=options.min_nvoxels)
    else:
        # Should never happen
        raise RuntimeError("No native ROI to get stats: %s" % str(roi))

def add_native_roi(rois, roi, name, threshold=0.5, log=sys.stdout):
    """
    Add an ROI in native (ASL) space
    """
    rois.append({"name" : name, "mask_native" : roi.data > threshold})
    log.write(" - %s...DONE\n" % name)

def add_struct_roi(rois, roi, name, ref, struct2asl, threshold=0.5, log=sys.stdout):
    """
    Add an ROI in structural space
    """
    log.write(" - %s..." % name)
    roi_native = _transform(roi, warp=None, ref=ref, premat=struct2asl)
    rois.append({"name" : name, "roi_struct" : roi, "roi_native" : roi_native, "mask_native" : roi_native.data > threshold})
    log.write("DONE\n")

def add_mni_roi(rois, roi, name, mni2struc, ref, struct2asl, threshold=0.5, log=sys.stdout):
    """
    Add an ROI in MNI space
    """
    log.write(" - %s..." % name)
    roi_native = _transform(roi, warp=mni2struc, ref=ref, postmat=struct2asl)
    rois.append({"name" : name, "roi_mni" : roi, "roi_native" : roi_native, "mask_native" : roi_native.data > threshold})
    log.write("DONE\n")

def apply_psf(array, psf):
    """
    Apply PSF blurring to an array (typically a binary or PV mask)

    The PSF is assumed to act only along the Z axis
    """
    if psf is None:
        return array

    # Make sure array is 4D
    array = array.astype(np.float)
    was_3d = False
    if array.ndim == 3:
        was_3d = True
        array = array[..., np.newaxis]

    # Detect Z dimension padding by comparing size of psf and data
    n_slices = psf.shape[0]
    padding_slices = n_slices - array.shape[2]
    if padding_slices < 0 or padding_slices % 2 != 0:
        raise ValueError("Invalid padding in psf: %i slices vs %i in data (difference must be even and > 0)" % (n_slices, array.shape[2]))
    padding_slices = int(padding_slices/2)
    array = np.pad(array, [(0, 0), (0, 0), (padding_slices, padding_slices), (0, 0)], 'edge')

    # Calculate mean along z direction for each (x, y, t) to demean volume
    zmean = np.expand_dims(np.mean(array, 2), 2)
    array = array - zmean

    # Apply blurring using multiplication in Fourier domain
    fftkern = fft(psf)[np.newaxis, np.newaxis, ..., np.newaxis]
    fftvol = fft(array, axis=2)

    # Get blurred volume in Image domain and add DC term back
    blurred_array = np.real(ifft(np.multiply(fftvol, fftkern), axis=2)) + zmean

    # Unpad and Unsqueeze extra dimension if original volume was 3D
    blurred_array = blurred_array[:, :, padding_slices:-padding_slices, :]
    if was_3d:
        blurred_array = np.squeeze(blurred_array, axis=3)
    return blurred_array

def add_struct_roi_set(rois, roi_set, names, ref, struct2asl, threshold=None, psf=None, log=sys.stdout):
    """
    Add an ROI set defined in structural space with optional PSF

    :param rois: Mapping from name to ROI array which will be updated
    :param roi_set: 4D Image where volumes define disjoint masks in MNI space
    :param names: Array of ROI names
    :param mni2struc: Warp image containing MNI->structural space warp
    :param struct2asl: Matrix for struct->ASL transformation
    :param threshol: Threshold for generating binary native space mask
    :param psf: Optional point spread function for binary ROIs. If specified, will be
                applied to ROIs and 'fuzzy' stats will be used
    """
    log.write(" - %s..." % ",".join(names))
    roi_set_native = _transform(roi_set, warp=None, ref=ref, premat=struct2asl)
    if threshold:
        mask_set_native = roi_set_native.data > threshold
    else:
        mask_set_native = roi_set_native.data
    fuzzy_set_native = apply_psf(mask_set_native, psf)
    rois.append({"names" : names, "roi_struct" : roi_set, "roi_native" : roi_set_native, "mask_native" : mask_set_native, "fuzzy_native" : fuzzy_set_native})
    log.write("DONE\n")

def add_mni_roi_set(rois, roi_set, names, mni2struc, ref, struct2asl, threshold=None, psf=None, log=sys.stdout):
    """
    Add an ROI set defined in MNI space with optional PSF

    :param rois: Mapping from name to ROI array which will be updated
    :param roi_set: 4D Image where volumes define disjoint masks in MNI space
    :param names: Array of ROI names
    :param mni2struc: Warp image containing MNI->structural space warp
    :param struct2asl: Matrix for struct->ASL transformation
    :param threshol: Threshold for generating binary native space mask
    :param psf: Optional point spread function for binary ROIs. If specified, will be
                applied to ROIs and 'fuzzy' stats will be used
    """
    log.write(" - %s..." % ",".join(names))
    roi_set_native = _transform(roi_set, warp=mni2struc, ref=ref, postmat=struct2asl)
    if threshold:
        mask_set_native = roi_set_native.data > threshold
    else:
        mask_set_native = roi_set_native.data
    fuzzy_set_native = apply_psf(mask_set_native, psf)
    rois.append({"names" : names, "roi_mni" : roi_set, "roi_native" : roi_set_native, "mask_native" : mask_set_native, "fuzzy_native" : fuzzy_set_native})
    log.write("DONE\n")

def add_roi_set_from_fsl_atlas(rois, mni2struc_warp, ref_img, struct2asl_mat, atlas_name, resolution=2, threshold=50, psf=None, log=sys.stdout):
    """
    Get ROIs from an FSL atlas

    :param rois: Mapping from name to ROI array which will be updated
    :param mni2struc_warp: Warp image containing MNI->structural space warp
    :param struct2asl_mat: Matrix for struct->ASL transformation
    :param atlas_name: Name of the FSL atlas
    :param resolution: Resolution in mm
    :param threshold: Threshold for probabilistic atlases
    :param psf: Optional point spread function for binary ROIs. If specified, will be
                applied to ROIs and 'fuzzy' stats will be used
    """
    log.write("\nAdding ROI set from standard atlas: %s (resolution=%imm, thresholding at %.2f)\n" % (atlas_name, resolution, threshold))
    registry = atlases.registry
    registry.rescanAtlases()
    desc = registry.getAtlasDescription(atlas_name)
    atlas = registry.loadAtlas(desc.atlasID, resolution=2)
    roi_set, names = [], []
    for label in desc.labels:
        roi_region = atlas.get(label=label)
        # When treating as an ROI set, convert to probability and do not threshold
        roi_set.append(roi_region.data / 100)
        names.append(label.name)
        if psf is None:
            add_mni_roi(rois, roi_region, label.name, mni2struc_warp, ref_img, struct2asl_mat, threshold=threshold)
    if psf is not None:
        roi_set = Image(np.stack(roi_set, axis=3), header=roi_region.header)
        add_mni_roi_set(rois, roi_set, names, mni2struc_warp, ref_img, struct2asl_mat, psf=psf)

def add_roi_set_from_mni_label_atlas(rois, mni2struc_warp, ref_img, struct2asl_mat, atlas_img, region_names, psf=None, log=sys.stdout):
    """
    Get ROIs from an atlas described by a label image in MNI space

    :param rois: Mapping from name to ROI array which will be updated
    :param mni2struc_warp: Warp image containing MNI->structural space warp
    :param struct2asl_mat: Matrix for struct->ASL transformation
    :param atlas_img: Atlas label image
    :param region_names: array of names for atlas regions
    :param psf: Optional point spread function for binary ROIs. If specified, will be
                applied to ROIs and 'fuzzy' stats will be used
    """
    log.write("\nAdding ROI set from MNI atlas image: %s\n" % (atlas_img.name))
    labels = [idx for idx in np.unique(atlas_img.data) if idx != 0]
    if len(labels) != len(region_names):
        region_names = ["Region %i" % label for label in labels]
    roi_set = []
    for name, label in zip(region_names, labels):
        roi = atlas_img.data.copy()
        roi_bin = (roi == label).astype(np.int)
        roi_set.append(roi_bin)
        roi_bin = Image(roi_bin, header=atlas_img.header)
        if psf is None:
            add_mni_roi(rois, roi_bin, name, mni2struc_warp, ref_img, struct2asl_mat, threshold=0.5)
    if psf is not None:
        roi_set = Image(np.stack(roi_set, axis=3), header=atlas_img.header)
        add_mni_roi_set(rois, roi_set, region_names, mni2struc_warp, ref_img, struct2asl_mat, psf=psf)

def get_perfusion_data(outdir, gm_pve_asl, wm_pve_asl, gm_thresh, wm_thresh, min_gm_thresh, min_wm_thresh, log=sys.stdout):
    brain_mask = Image(os.path.join(outdir, "mask"))
    perfusion_data = [
        {
            "suffix" : "", 
            "f" : Image(os.path.join(outdir, "perfusion_calib")), 
            "var" : Image(os.path.join(outdir, "perfusion_var_calib")),
            "mask" : brain_mask.data,
        },
    ]
    if os.path.isdir(os.path.join(outdir, "pvcorr")):
        log.write(" - Found partial volume corrected results - will mask ROIs using 'base' GM/WM masks (PVE thresholds: %.2f / %.2f)\n" % (min_gm_thresh, min_wm_thresh))
        perfusion_data.extend([
            {
                "suffix" : "_gm", 
                "f" : Image(os.path.join(outdir, "pvcorr", "perfusion_calib")), 
                "var" : Image(os.path.join(outdir, "pvcorr", "perfusion_var_calib")),
                "mask" : np.logical_and(brain_mask.data, gm_pve_asl.data > min_gm_thresh),
            },
            {
                "suffix" : "_wm", 
                "f" : Image(os.path.join(outdir, "pvcorr", "perfusion_wm_calib")), 
                "var" : Image(os.path.join(outdir, "pvcorr", "perfusion_wm_var_calib")),
                "mask" : np.logical_and(brain_mask.data, wm_pve_asl.data > min_wm_thresh),
            },
        ])
    else:
        log.write(" - No partial volume corrected results - will mask ROIs using 'pure' GM/WM masks (PVE thresholds: %.2f / %.2f)\n" % (gm_thresh, wm_thresh))
        perfusion_data.extend([
            {
                "suffix" : "_gm",
                "f" : Image(os.path.join(outdir, "perfusion_calib")), 
                "var" : Image(os.path.join(outdir, "perfusion_var_calib")),
                "mask" : np.logical_and(brain_mask.data, gm_pve_asl.data > gm_thresh),
            },
            {
                "suffix" : "_wm",
                "f" : Image(os.path.join(outdir, "perfusion_calib")), 
                "var" : Image(os.path.join(outdir, "perfusion_var_calib")),
                "mask" : np.logical_and(brain_mask.data, wm_pve_asl.data > wm_thresh),
            },
        ])
    return perfusion_data

def get_arrival_data(outdir, gm_pve_asl, wm_pve_asl, gm_thresh, wm_thresh, min_gm_thresh, min_wm_thresh, log=sys.stdout):
    brain_mask = Image(os.path.join(outdir, "mask"))
    f = Image(os.path.join(outdir, "perfusion_calib"))
    f_var = Image(os.path.join(outdir, "perfusion_var_calib"))
    f_good = np.logical_and(np.isfinite(f.data), f_var.data > 0)

    # Note that for the arrival mask we also remove voxels which will be
    # eliminated from the perfusion stats because of nan/inf values or zero
    # variances. This is a bit of a hack but should help ensure that the
    # voxel set is consistent between the two measures. Perfusion can have
    # invalid values in voxels where the arrival time is valid because of
    # voxelwise calibration, however we do not expect the reverse to occur
    arrival_data = [
        {
            "suffix" : "_arrival",
            "f" : Image(os.path.join(outdir, "arrival")),
            "var" : Image(os.path.join(outdir, "arrival_var")),
            "mask" : np.logical_and(brain_mask.data, f_good)
        },
    ]
    if os.path.isdir(os.path.join(outdir, "pvcorr")):
        arrival_data.extend([
            {
                "suffix" : "_arrival_gm",
                "f" : Image(os.path.join(outdir, "pvcorr", "arrival")),
                "var" : Image(os.path.join(outdir, "pvcorr", "arrival_var")),
                "mask" : np.logical_and(brain_mask.data, gm_pve_asl.data > min_gm_thresh),
            },
            {
                "suffix" : "_arrival_wm",
                "f" : Image(os.path.join(outdir, "pvcorr", "arrival_wm")),
                "var" : Image(os.path.join(outdir, "pvcorr", "arrival_wm_var")),
                "mask" : np.logical_and(brain_mask.data, wm_pve_asl.data > min_wm_thresh),
            },
        ])
    else:
        arrival_data.extend([
            {
                "suffix" : "_arrival_gm",
                "f" : Image(os.path.join(outdir, "arrival")),
                "var" : Image(os.path.join(outdir, "arrival_var")),
                "mask" : np.logical_and(
                    np.logical_and(brain_mask.data, gm_pve_asl.data > gm_thresh),
                    f_good
                )
            },
            {
                "suffix" : "_arrival_wm",
                "f" : Image(os.path.join(outdir, "arrival")),
                "var" : Image(os.path.join(outdir, "arrival_var")),
                "mask" : np.logical_and(
                    np.logical_and(brain_mask.data, wm_pve_asl.data > wm_thresh),
                    f_good
                )
            },
        ])
    return arrival_data

def main():
    options = ArgumentParser().parse_args()
    if options.oxasl_output is None:
        sys.stderr.write("oxford_asl output directory must be specified")
        sys.exit(1)

    print("Regionwise analysis\n")
    print(" - Using oxford_asl output in %s" % options.oxasl_output)

    # Get reference and transformation data from oxford_asl and fsl_anat output
    outdir = options.oxasl_output
    asl_ref = Image(os.path.join(outdir, "perfusion"))

    gm_pve, wm_pve, csf_pve = None, None, None
    if options.fslanat is not None:
        print(" - Using fsl_anat output in %s" % options.fslanat)
        struc_ref = Image(os.path.join(options.fslanat, "T1"))
        gm_pve = Image(os.path.join(options.fslanat, "T1_fast_pve_1"))
        wm_pve = Image(os.path.join(options.fslanat, "T1_fast_pve_2"))
        csf_pve = Image(os.path.join(options.fslanat, "T1_fast_pve_0"))
        struct2mni_warp = Image(os.path.join(options.fslanat, "T1_to_MNI_nonlin_coeff"))
    elif options.struc is not None and (options.struc2std is not None or options.std2struc is not None):
        print(" - Using manually specified structural data")
        struc_ref = Image(options.struc)
        if options.std2struc is not None:
            print(" - std->struc transformation warp provided directly")
            mni2struc_warp = Image(options.std2struc)
        elif options.struc2std is not None:
            struct2mni_warp = Image(options.struc2std)
        if options.gm_pve is not None:
            gm_pve = Image(options.gm_pve)
        if options.wm_pve is not None:
            wm_pve = Image(options.wm_pve)
        if options.csf_pve is not None:
            csf_pve = Image(options.csf_pve)
    else:
        sys.stderr.write("Either --fslanat must be specified or all of --struc, --gm-pve, --wm-pve --csf-pve and --struc2std/--std2struc \n")
        sys.exit(1)

    if options.fslanat is not None or options.std2struc is None:
        print(" - Generating std->struc transformation by inverting struc->std transformation warp")
        mni2struc_warp = fsl.invwarp(struct2mni_warp, struc_ref, out=fsl.LOAD)["out"]

    asl2struc_filename = options.asl2struc
    if asl2struc_filename is None:
        asl2struc_filename = os.path.join(outdir, "asl2struct.mat")
    with open(asl2struc_filename) as asl2struct_file:
        print(" - Loading ASL->struc transformation matrix from %s" % asl2struc_filename)
        asl2struct_mat = np.array([[float(v) for v in line.split()] for line in asl2struct_file.readlines()])
        struct2asl_mat = np.linalg.inv(asl2struct_mat)

    if options.psf:
        sys.stdout.write(f" - Loading PSF {options.psf} for ROI sets: ")
        options.psf = np.loadtxt(options.psf)
        print("DONE - %i Z slice values in PSF" % len(options.psf))

    # Look for PVC or non-PVC results
    print("\nLoading perfusion images")
    if options.native_pves:
        print("Using native space PVEs.")
        gm_pve_asl = gm_pve
        wm_pve_asl = wm_pve
    else:
        print("Transforming PVEs from structural to native space.")
        gm_pve_asl = _transform(gm_pve, warp=None, ref=asl_ref, premat=struct2asl_mat)
        wm_pve_asl = _transform(wm_pve, warp=None, ref=asl_ref, premat=struct2asl_mat)
    perfusion_data = get_perfusion_data(outdir, gm_pve_asl, wm_pve_asl, options.gm_thresh, options.wm_thresh, options.min_gm_thresh, options.min_wm_thresh)
    if options.add_arrival:
        print(" - Also generating stats for arrival data")
        perfusion_data += get_arrival_data(outdir, gm_pve_asl, wm_pve_asl, options.gm_thresh, options.wm_thresh, options.min_gm_thresh, options.min_wm_thresh)

    rois = []
    print("\nLoading generic ROIs")
    if gm_pve is not None and not options.native_pves:
        add_struct_roi(rois, gm_pve, "%i%%+GM" % (options.min_gm_thresh*100), ref=asl_ref, struct2asl=struct2asl_mat, threshold=options.min_gm_thresh)
        add_struct_roi(rois, gm_pve, "%i%%+GM" % (options.gm_thresh*100), ref=asl_ref, struct2asl=struct2asl_mat, threshold=options.gm_thresh)
    elif gm_pve is not None:
        add_native_roi(rois, gm_pve_asl, "%i%%+GM" % (options.min_gm_thresh*100), threshold=options.min_gm_thresh)
        add_native_roi(rois, gm_pve_asl, "%i%%+GM" % (options.gm_thresh*100), threshold=options.gm_thresh)
    if wm_pve is not None and not options.native_pves:
        add_struct_roi(rois, wm_pve, "%i%%+WM" % (options.min_wm_thresh*100), ref=asl_ref, struct2asl=struct2asl_mat, threshold=options.min_wm_thresh)
        add_struct_roi(rois, wm_pve, "%i%%+WM" % (options.wm_thresh*100), ref=asl_ref, struct2asl=struct2asl_mat, threshold=options.wm_thresh)
    elif wm_pve is not None:
        add_native_roi(rois, wm_pve_asl, "%i%%+WM" % (options.min_wm_thresh*100), threshold=options.min_wm_thresh)
        add_native_roi(rois, wm_pve_asl, "%i%%+WM" % (options.wm_thresh*100), threshold=options.wm_thresh)

    if gm_pve is not None and wm_pve is not None and csf_pve is not None and not options.native_pves:
        print("\nLoading tissue PV ROI set")
        roi_set = Image(np.stack([gm_pve.data, wm_pve.data, csf_pve.data], axis=-1), header=gm_pve.header)
        add_struct_roi_set(rois, roi_set, ["GM PV", "WM PV", "CSF PV"], ref=asl_ref, struct2asl=struct2asl_mat, psf=options.psf)

    # Add ROIs from command line
    print("\nLoading user-specified ROIs")
    for fname in options.roi_native:
        add_native_roi(rois, Image(fname), os.path.basename(fname).split(".")[0])
    for fname in options.roi_struct:
        add_struct_roi(rois, Image(fname), os.path.basename(fname).split(".")[0], asl_ref, struct2asl_mat)

    print("\nLoading user-specified label atlases")
    for idx, fname in enumerate(options.add_mni_atlas):
        if idx < len(options.add_mni_atlas_labels):
            with open(options.add_mni_atlas_labels[idx]) as f:
                names = [l.strip() for l in f.readlines()]
        else:
            names = [os.path.basename(fname).split(".")[0],]
        add_roi_set_from_mni_label_atlas(rois, mni2struc_warp, asl_ref, struct2asl_mat, Image(fname), names, psf=options.psf)

    # Add ROIs from standard atlases
    if options.add_standard_atlases:
        # FIXME combine as single ROI set?
        add_roi_set_from_fsl_atlas(rois, mni2struc_warp, asl_ref, struct2asl_mat, "harvardoxford-cortical", psf=options.psf)
        add_roi_set_from_fsl_atlas(rois, mni2struc_warp, asl_ref, struct2asl_mat, "harvardoxford-subcortical", psf=options.psf)

    # Get stats in each ROI. Add name to stats dict to make TSV output easier
    print("\nGetting stats - minimum of %i voxels to report in region" % options.min_nvoxels)
    for item in perfusion_data:
        suffix = item["suffix"]
        os.makedirs(options.output, exist_ok=True)
        writer = None
        # Save TSV stats output. Note we give TSV file a CSV extension to make Excel happier
        with open(os.path.join(options.output, "%s%s.csv" % (options.output_prefix, suffix)), mode="w", newline='') as tsv_file:
            for roi in rois:
                roi["stats"] = {"name" : roi.get("name", None), "names" : roi.get("names", None)}
                get_stats(roi, item, options)

                if writer is None:
                    writer = csv.DictWriter(tsv_file, fieldnames=list(roi["stats"].keys()))
                    writer.writeheader()

                # The ROI might just be a single ROI or it might be a 'fuzzy set' - need
                # to handle both cases
                roi_stats = dict(roi["stats"])
                if roi_stats["name"] is not None:
                    all_roi_stats = [roi_stats]
                else:
                    all_roi_stats = []
                    for idx in range(len(roi_stats["names"])):
                        roi_stats_copy = dict(roi_stats)
                        for k in list(roi_stats.keys()):
                            if roi_stats_copy[k] is not None:
                                roi_stats_copy[k] = roi_stats_copy[k][idx]
                        roi_stats_copy["name"] = roi_stats_copy.pop("names")
                        all_roi_stats.append(roi_stats_copy)

                # Round floats to 2 d.p. - should be fine for calibrated perfusion values
                for roi_stats in all_roi_stats:
                    rounded = dict(roi_stats)
                    for k, v in roi_stats.items():
                        if isinstance(v, (float, np.float, np.float32, np.float64)):
                            rounded[k] = "%.2f" % v
                    writer.writerow(rounded)

    # Save output masks/PVE maps
    for roi in rois:
        if "name" in roi:
            names = [roi["name"]]
        else:
            names = roi["names"]
        for idx, name in enumerate(names):
            fname = name.replace(" ", "_").replace(",", "").lower() + ".nii.gz"
            if options.save_native_rois:
                if "roi_native" in roi:
                    _write_nii(roi["roi_native"], os.path.join(options.output, "rois_native", fname), vol=idx)
                if "mask_native" in roi:
                    _write_nii(roi["mask_native"], os.path.join(options.output, "masks_native", fname), header=asl_ref.header, vol=idx)
                if "fuzzy_native" in roi:
                    _write_nii(roi["fuzzy_native"], os.path.join(options.output, "fuzzy_native", fname), header=asl_ref.header, vol=idx)
            if options.save_struct_rois and "roi_struct" in roi:
                _write_nii(roi["roi_struct"], os.path.join(options.output, "rois_struct", fname), vol=idx)
            if options.save_mni_rois and "roi_mni" in roi:
                _write_nii(roi["roi_mni"], os.path.join(options.output, "rois_mni", fname), vol=idx)

    print("\nDONE - Output in %s" % options.output)

if __name__ == "__main__":
    main()
