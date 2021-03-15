#!/bin/env python
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

class ArgumentParser(argparse.ArgumentParser):
    """
    ArgumentParser for program options
    """

    def __init__(self, **kwargs):
        argparse.ArgumentParser.__init__(self, prog="oxford_asl_roi_stats", add_help=False, **kwargs)
        self.add_argument("--oxasl-output", required=True,
                          help="OXFORD_ASL or OXASL output directory")
        self.add_argument("--fslanat",
                          help="FSL_ANAT output directory - if not specified --struc --gm-pve --wm-pve and --struc2std must be given")
        self.add_argument("--struc", "-s", help="Structural space reference image - ignored if --fslanat is given")
        self.add_argument("--gm-pve", help="GM PVE in structural space - ignored if --fslanat is given")
        self.add_argument("--wm-pve", help="WM PVE in structural space - ignored if --fslanat is given")
        self.add_argument("--struc2std", help="Structural -> standard space nonlinear warp map - ignored if --fslanat is given")
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
        self.add_argument("--save-native-rois", action="store_true", default=False,
                          help="Save ROIs in native (ASL) space")
        self.add_argument("--save-native-masks", action="store_true", default=False,
                          help="Save binary masks in native (ASL) space")

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
        ret = Image((ret.data > output_roi_threshold).astype(np.int), header=ret.header)
    return ret

def _write_nii(img, fname, header=None):
    os.makedirs(os.path.dirname(fname), exist_ok=True)
    if isinstance(img, Image):
        img.save(fname)
    else:
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

def i2(val, var):
    """ I^2 Measure of heterogenaity """
    prec = 1 / var
    prec[~np.isfinite(prec)] = 0
    n = len(val)
    mu_bar = mean_invvarweighted(val, var)

    #out = []
    Q = np.sum(prec * (val - mu_bar)**2)
    #H = np.sqrt(Q/(n - 1))
    if Q == 0:
        i2 = 0
    else:
        i2 = (Q-(n-1))/Q

    # Negative values map to 0 (see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC192859/)
    i2 = max(i2, 0)
    #tau2_DL = (Q - n + 1)/(sum(w) - sum([x**2 for x in w])/sum(w))
    #tau2_DL = max(0, tau2_DL)
    #out.extend([mu_bar,Q,H,I2,tau2_DL])

    #out = pd.DataFrame(out)
    #out['index']=['Weighted_Mean', 'Q', 'H', 'I2', 'tau2_DL']
    #out = out.set_index('index', drop=True)

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

def get_stats(stats, img, var_img, roi, suffix="", ignore_nan=True, ignore_inf=True, min_nvoxels=10, mask=None):
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

    effective_roi = roi
    if ignore_nan:
        effective_roi = np.logical_and(effective_roi, ~np.isnan(img))
    if ignore_inf:
        effective_roi = np.logical_and(effective_roi, np.isfinite(img))
    if mask is not None:
        effective_roi = np.logical_and(effective_roi, mask)

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

def add_native_roi(rois, roi, name, log=sys.stdout):
    """
    Add an ROI in native (ASL) space
    """
    rois.append({"name" : name, "mask_native" : roi.data})
    log.write(" - %s...DONE" % name)

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
    roi_native = _transform(roi, mni2struc, ref=ref, postmat=struct2asl)
    rois.append({"name" : name, "roi_mni" : roi, "roi_native" : roi_native, "mask_native" : roi_native.data > threshold})
    log.write("DONE\n")
 
def add_rois_from_fsl_atlas(rois, mni2struc_warp, ref_img, struct2asl_mat, atlas_name, resolution=2, threshold=50, log=sys.stdout):
    """
    Get ROIs from an FSL atlas
    
    :param rois: Mapping from name to ROI array which will be updated
    :param mni2struc_warp: Warp image containing MNI->structural space warp
    :param struct2asl_mat: Matrix for struct->ASL transformation
    :param atlas_name: Name of the FSL atlas
    :param resolution: Resolution in mm
    :param threshold: Threshold for probabilistic atlases
    """
    log.write("\nAdding ROIs from standard atlas: %s (resolution=%imm, thresholding at %.2f)\n" % (atlas_name, resolution, threshold))
    registry = atlases.registry
    registry.rescanAtlases()
    desc = registry.getAtlasDescription(atlas_name)
    atlas = registry.loadAtlas(desc.atlasID, resolution=2)
    for label in desc.labels:
        add_mni_roi(rois, atlas.get(label=label), label.name, mni2struc_warp, ref_img, struct2asl_mat, threshold=threshold)

def add_rois_from_mni_atlas(rois, mni2struc_warp, ref_img, struct2asl_mat, atlas_img, region_names, log=sys.stdout):
    """
    Get ROIs from an atlas described by a label image in MNI space
 
    :param rois: Mapping from name to ROI array which will be updated
    :param mni2struc_warp: Warp image containing MNI->structural space warp
    :param struct2asl_mat: Matrix for struct->ASL transformation
    :param atlas_img: Atlas label image
    :param region_names: Atlas label image
    :param resolution: Resolution in mm
    """
    log.write("\nAdding ROIs from MNI atlas image: %s\n" % (atlas_img.name))
    labels = [idx for idx in np.unique(atlas_img.data) if idx != 0]
    if len(labels) != len(region_names):
        region_names = ["Region %i" % label for label in labels]
    for idx, label in enumerate(labels):
        roi = atlas_img.data.copy()
        roi[roi != label] = 0
        name = region_names[idx]
        add_mni_roi(rois, Image(roi, header=atlas_img.header), name, mni2struc_warp, ref_img, struct2asl_mat, threshold=0.5)

def get_perfusion_data(outdir, gm_pve_asl, wm_pve_asl, gm_thresh, wm_thresh, min_gm_thresh, min_wm_thresh, log=sys.stdout):
    brain_mask_prepvc = Image(os.path.join(outdir, "mask"))
    perfusion_data = [
        {
            "suffix" : "", 
            "f" : Image(os.path.join(outdir, "perfusion_calib")), 
            "var" : Image(os.path.join(outdir, "perfusion_var_calib")),
            "mask" : brain_mask_prepvc.data,
        },
    ]
    if os.path.isdir(os.path.join(outdir, "pvcorr")):
        log.write(" - Found partial volume corrected results - will mask ROIs using 'base' GM/WM masks (PVE thresholds: %.2f / %.2f)\n" % (min_gm_thresh, min_wm_thresh))
        brain_mask_pvc = Image(os.path.join(outdir, "mask_pvcorr"))
        perfusion_data.extend([
            {
                "suffix" : "_gm", 
                "f" : Image(os.path.join(outdir, "pvcorr", "perfusion_calib")), 
                "var" : Image(os.path.join(outdir, "pvcorr", "perfusion_var_calib")),
                "mask" : np.logical_and(brain_mask_pvc.data, gm_pve_asl.data > min_gm_thresh),
            },
            {
                "suffix" : "_wm", 
                "f" : Image(os.path.join(outdir, "pvcorr", "perfusion_wm_calib")), 
                "var" : Image(os.path.join(outdir, "pvcorr", "perfusion_wm_var_calib")),
                "mask" : np.logical_and(brain_mask_pvc.data, wm_pve_asl.data > min_wm_thresh),
            },
        ])
    else:
        log.write(" - No partial volume corrected results - will mask ROIs using 'pure' GM/WM masks (PVE thresholds: %.2f / %.2f)\n" % (gm_thresh, wm_thresh))
        perfusion_data.extend([
            {
                "suffix" : "_gm",
                "f" : Image(os.path.join(outdir, "perfusion_calib")), 
                "var" : Image(os.path.join(outdir, "perfusion_var_calib")),
                "mask" : np.logical_and(brain_mask_prepvc.data, gm_pve_asl.data > gm_thresh),
            },
            {
                "suffix" : "_wm",
                "f" : Image(os.path.join(outdir, "perfusion_calib")), 
                "var" : Image(os.path.join(outdir, "perfusion_var_calib")),
                "mask" : np.logical_and(brain_mask_prepvc.data, wm_pve_asl.data > wm_thresh),
            },
        ])
    return perfusion_data

def main():
    options = ArgumentParser().parse_args()
    if options.oxasl_output is None:
        sys.stderr.write("oxford_asl output directory must be specified")
        sys.exit(1)

    print("Regionwise analysis\n")
    print(" - Using oxford_asl output in %s" % options.oxasl_output)
    print(" - Using fsl_anat output in %s" % options.fslanat)

    # Get reference and transformation data from oxford_asl and fsl_anat output
    outdir = os.path.join(options.oxasl_output, "native_space")
    asl_ref = Image(os.path.join(outdir, "perfusion"))

    if options.fslanat is not None:
        struc_ref = Image(os.path.join(options.fslanat, "T1"))
        gm_pve = Image(os.path.join(options.fslanat, "T1_fast_pve_1"))
        wm_pve = Image(os.path.join(options.fslanat, "T1_fast_pve_2"))
        struct2mni_warp = Image(os.path.join(options.fslanat, "T1_to_MNI_nonlin_coeff"))
    elif options.struc is not None and options.gm_pve is not None and options.wm_pve is not None and options.struc2std is not None:
        struc_ref = Image(options.struc)
        gm_pve = Image(options.gm_pve)
        wm_pve = Image(options.wm_pve)
        struct2mni_warp = Image(options.struc2std)
    else:
        sys.stderr.write("Either --fslanat must be specified or all of --struc, --gm-pve, --wm-pve and --struc2std")
        sys.exit(1)

    mni2struc_warp = fsl.invwarp(struct2mni_warp, struc_ref, out=fsl.LOAD)["out"]
    with open(os.path.join(outdir, "asl2struct.mat")) as asl2struct_file:
        asl2struct_mat = np.array([[float(v) for v in line.split()] for line in asl2struct_file.readlines()])
        struct2asl_mat = np.linalg.inv(asl2struct_mat)

    # Look for PVC or non-PVC results
    print("\nLoading perfusion images")
    gm_pve_asl = _transform(gm_pve, warp=None, ref=asl_ref, premat=struct2asl_mat)
    wm_pve_asl = _transform(wm_pve, warp=None, ref=asl_ref, premat=struct2asl_mat)
    perfusion_data = get_perfusion_data(outdir, gm_pve_asl, wm_pve_asl, options.gm_thresh, options.wm_thresh, options.min_gm_thresh, options.min_wm_thresh)

    rois = []
    print("\nLoading generic ROIs")
    add_struct_roi(rois, gm_pve, "%i%%+GM" % (options.min_gm_thresh*100), ref=asl_ref, struct2asl=struct2asl_mat, threshold=options.min_gm_thresh)
    add_struct_roi(rois, wm_pve, "%i%%+WM" % (options.min_wm_thresh*100), ref=asl_ref, struct2asl=struct2asl_mat, threshold=options.min_wm_thresh)
    add_struct_roi(rois, gm_pve, "%i%%+GM" % (options.gm_thresh*100), ref=asl_ref, struct2asl=struct2asl_mat, threshold=options.gm_thresh)
    add_struct_roi(rois, wm_pve, "%i%%+WM" % (options.wm_thresh*100), ref=asl_ref, struct2asl=struct2asl_mat, threshold=options.wm_thresh)

    # Add ROIs from command line
    print("\nLoading user-specified ROIs")
    for fname in options.roi_native:
        add_native_roi(rois, Image(fname), os.path.basename(fname).split(".")[0])
    for fname in options.roi_struct:
        add_struct_roi(rois, Image(fname), os.path.basename(fname).split(".")[0], asl_ref, struct2asl_mat)
    for idx, fname in enumerate(options.add_mni_atlas):
        if idx < len(options.add_mni_atlas_labels):
            with open(options.add_mni_atlas_labels[idx]) as f:
                names = [l.strip() for l in f.readlines()]
        else:
            names = [os.path.basename(fname).split(".")[0],]
        add_rois_from_mni_atlas(rois, mni2struc_warp, asl_ref, struct2asl_mat, Image(fname), names)

    # Add ROIs from standard atlases
    if options.add_standard_atlases:
        add_rois_from_fsl_atlas(rois, mni2struc_warp, asl_ref, struct2asl_mat, "harvardoxford-cortical")
        add_rois_from_fsl_atlas(rois, mni2struc_warp, asl_ref, struct2asl_mat, "harvardoxford-subcortical")

    # Get stats in each ROI. Add name to stats dict to make TSV output easier
    print("\nGetting stats - minimum of %i voxels to report in region" % options.min_nvoxels)
    for item in perfusion_data:
        suffix = item["suffix"]
        os.makedirs(options.output, exist_ok=True)
        writer = None
        # Save TSV stats output. Note we give TSV file a CSV extension to make Excel happier
        with open(os.path.join(options.output, "region_analysis%s.csv" % suffix), mode="w", newline='') as tsv_file:
            for roi in rois:
                roi["stats"] = {"name" : roi["name"]}
                get_stats(roi["stats"], item["f"].data, item["var"].data, roi["mask_native"], mask=item["mask"], min_nvoxels=options.min_nvoxels)
            
                if writer is None:
                    writer = csv.DictWriter(tsv_file, fieldnames=list(roi["stats"].keys()))
                    writer.writeheader()

                # Round floats to 2 d.p. - should be fine for calibrated perfusion values
                rounded = dict(roi["stats"])
                for k, v in roi["stats"].items():
                    if isinstance(v, (float, np.float, np.float32, np.float64)):
                        rounded[k] = "%.2f" % v
                writer.writerow(rounded)

    # Save output masks/PVE maps
    for roi in rois:
        fname = roi["name"].replace(" ", "_").replace(",", "").lower() + ".nii.gz"
        if options.save_native_rois and "roi_native" in roi:
            _write_nii(roi["roi_native"], os.path.join(options.output, "rois_native", fname))
        if options.save_native_masks and "mask_native" in roi:
            _write_nii(roi["mask_native"], os.path.join(options.output, "masks_native", fname), header=asl_ref.header)
        if options.save_mni_rois and "roi_mni" in roi:
            _write_nii(roi["roi_mni"], os.path.join(options.output, "rois_mni", fname))

    print("\nDONE - Output in %s" % options.output)

if __name__ == "__main__":
    main()
