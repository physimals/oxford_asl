#!/usr/bin/env python
"""
Basic classes for ASL data processing
"""

import os
import sys
import shutil
from optparse import OptionParser, OptionGroup

import numpy as np
import nibabel as nib

import asl.fslwrap as fsl

class AslOptionParser(OptionParser):
    """
    OptionParser which incorporates standard ASL options
    """
    def __init__(self, asl_fname_opt="-i", ignore_stdopt=None, *args, **kwargs):
        """
        :param asl_fname_opt: Option which identifies the filename of ASL data
        :param ignore_stdopt: List of standard options to ignore (because they are not useful 
                              in this context)
        """
        OptionParser.__init__(self, *args, **kwargs)
        self.ignore_stdopt = []
        if ignore_stdopt is not None:
            self.ignore_stdopt = ignore_stdopt

        g = OptionGroup(self, "Input data")
        self._add(g, asl_fname_opt, dest="asldata", help="ASL data file")
        self._add(g, "--order", dest="order", help="Data order as sequence of 2 or 3 characters: t=TIs/PLDs, r=repeats, p/P=TC/CT pairs. First character is fastest varying", default="prt")
        self._add(g, "--tis", dest="tis", help="TIs as comma-separated list")
        self._add(g, "--plds", dest="plds", help="PLDs as comma-separated list")
        self._add(g, "--nrpts", dest="nrpts", help="Fixed number of repeats per TI", default=None)
        self._add(g, "--rpts", dest="rpts", help="Variable repeats as comma-separated list, one per TI", default=None)
        self._add(g, "--t1", dest="t1", help="T1 image")
        self._add(g, "--t1b", dest="t1b", help="Blood t1", type=float, default=1.65)
        self._add(g, "--te", dest="te", help="Echo time", type=float, default=0.012)
        self._add(g, "--alpha", dest="alpha", help="", type=float, default=0.98)
        self._add(g, "--lambda", dest="lam", help="", type=float, default=0.9)
        self._add(g, "--debug", dest="debug", help="Debug mode", action="store_true", default=False)
        self.add_option_group(g)

        g = OptionGroup(self, "Preprocessing")
        self._add(g, "--diff", dest="diff", help="Perform tag-control subtraction", action="store_true", default=False)
        self._add(g, "--smooth", dest="smooth", help="Spatially smooth data", action="store_true", default=False)
        self._add(g, "--fwhm", dest="fwhm", help="FWHM for spatial filter kernel", type="float", default=6)
        self._add(g, "--mc", dest="mc", help="Motion correct data", action="store_true", default=False)
        self._add(g, "--reorder", dest="reorder", help="Re-order data in specified order")
        self.add_option_group(g)

    def _add(self, g, name, *args, **kwargs):
        if name not in self.ignore_stdopt:
            g.add_option(name, *args, **kwargs)

class AslImage(fsl.Image):
  
    def __init__(self, name, order="prt", ntis=None, tis=None, plds=None, nrpts=None, rpts=None, **kwargs):
        fsl.Image.__init__(self, name, **kwargs)
        if self.ndim != 4:
            raise RuntimeError("4D data expected")
        self.nv = self.shape[3]

        self.order = order
        if "p" in order or "P" in order:
            self.ntc = 2
            self.tagfirst = "p" in self.order
        else:
            self.ntc = 1
            self.tagfirst = False

        if tis is not None and plds is not None:
            raise RuntimeError("Cannot specify PLDs and TIs at the same time")
        elif plds is not None:
            tis = plds
            self.plds = True
        else:
            self.plds = False

        if ntis is None and tis is None:
            raise RuntimeError("Number of TIs not specified")
        elif tis is not None:
            if isinstance(tis, basestring): tis = [float(ti) for ti in tis.split(",")]
            ntis = len(tis)
            if ntis is not None and len(tis) != ntis:
                raise RuntimeError("Number of TIs: %i, but list of %i TIs given" % (ntis, len(tis)))
        self.tis = tis
        self.ntis = int(ntis)
        
        if nrpts is not None and rpts is not None:
            raise RuntimeError("Cannot specify both fixed and variable numbers of repeats")        
        elif nrpts is None and rpts is None:
            # Calculate fixed number of repeats 
            if self.nv % (self.ntc * self.ntis) != 0:
                raise RuntimeError("Data contains %i volumes, inconsistent with %i TIs and TC pairs" % (self.nv, self.ntis))        
            rpts = [self.nv / (self.ntc * self.ntis)] * self.ntis
        elif nrpts is not None:
            nrpts = int(nrpts)
            if nrpts * self.ntis * self.ntc != self.nv:
                raise RuntimeError("Data contains %i volumes, inconsistent with %i TIs and %i repeats" % (self.nv, self.ntis, nrpts))
            rpts = [nrpts] * self.ntis
        else:
            if isinstance(rpts, basestring): rpts = [int(rpt) for rpt in rpts.split(",")]
            if len(rpts) != self.ntis:
                raise RuntimeError("%i TIs, but only %i variable repeats" % (self.ntis, len(rpts)))        
            elif sum(rpts) * self.ntc != self.nv:
                raise RuntimeError("Data contains %i volumes, inconsistent with total number of variable repeats" % self.nv)        

        self.rpts = rpts
        
    def _get_order_idx(self, order, tag, ti, rpt):
        idx = 0
        first = True
        for comp in order[::-1]:
            #print("comp: %s" % comp)
            if not first:
                idx *= self._get_ncomp(comp, ti)
                #print("Multiplied by %i" % self._get_ncomp(comp, ti))
            idx += self._get_comp(comp, tag, ti, rpt)
            #print("Added %i" % self._get_comp(comp, tag, ti, rpt))
            first = False
        return idx

    def _get_comp(self, comp_id, tag, ti, rpt):
        ret = {"t": ti, "r" : rpt, "p" : tag, "P" : 1-tag}
        if comp_id in ret: 
            return ret[comp_id]
        else:
            raise RuntimeError("Unknown ordering character: %s" % comp_id)

    def _get_ncomp(self, comp_id, ti):
        ret = {"t": self.ntis, "r" : self.rpts[ti], "p" : 2, "P" : 2}
        if comp_id in ret: 
            return ret[comp_id]
        else:
            raise RuntimeError("Unknown ordering character: %s" % comp_id)

    def reorder(self, out_order):
        """
        Re-order ASL data 

        The order is defined by a string in which
        r=repeats, p=tag-control pairs, P=control-tag pairs and t=tis/plds.
        The first character is the fastest varying

        So for a data set with 3 TIs and 2 repeats an order of "ptr" would be:
        TC (TI1), TC (TI2), TC (TI3), TC(TI1, repeat 2), TC(TI2 repeat 2), etc.
        """
        if self.ntc == 1 and ("p" in out_order or "P" in out_order):
            raise RuntimeError("Data contains TC pairs but output order does not")
        elif self.ntc == 2 and ("p" not in out_order and "P" not in out_order):
            raise RuntimeError("Output order contains TC pairs but input data  does not")

        #print("reordering from %s to %s" % (self.order, out_order))
        output_data = np.zeros(self.shape)
        input_data = self.data()
        tags = range(self.ntc)
        for ti in range(self.ntis):
            for rpt in range(self.rpts[ti]):
                for tag in tags:
                    #print("ti=%i, rpt=%i, tag=%i" % (ti, rpt, tag))
                    in_idx = self._get_order_idx(self.order, tag, ti, rpt)
                    #print("Input (%s) index %i" % (self.order, in_idx))
                    out_idx = self._get_order_idx(out_order, tag, ti, rpt)
                    #print("Output (%s) index %i" % (out_order, out_idx))
                    output_data[:,:,:,out_idx] = input_data[:,:,:,in_idx]
                    #print("")
        return AslImage(self.ipath + "_reorder", data=output_data,
                        order=out_order, tis=self.tis, ntis=self.ntis, rpts=self.rpts,
                        base=self)

    def single_ti(self, ti_idx, order=None):
            
        if order is None:
            if self.ntc == 2: order = "pr"
            else: order = "r"
        elif "t" in order:
            order = order.remove("t")
        order = order + "t"

        # Re-order so that TIs are together
        reordered = self.reorder(order)

        # Find the start index for this TI and the number of times it was repeated
        start = 0
        for idx in range(ti_idx):
            start += self.rpts[idx]*self.ntc
        nrpts = self.rpts[ti_idx]
        nvols = nrpts * self.ntc
        output_data = reordered.data()[:,:,:,start:start+nvols]
        if self.tis is not None:
            tis = [self.tis[ti_idx],]
        else:
            tis = None
        return AslImage(self.ipath + "_ti%i" % ti_idx, data=output_data,
                        order=order, tis=tis, ntis=1, nrpts=nrpts, base=self)

    def diff(self):
        """
        Perform tag-control differencing. Data is assumed to be ordered so that
        tag-control pairs are together
        """
        if "p" not in self.order and "P" not in self.order:
            # Already differenced
            output_data = self.data()
        else:
            output_data = np.zeros(list(self.shape[:3]) + [self.nv/2])

            # Re-order so that TC pairs are together with the tag first
            out_order = self.order.replace("p", "").replace("P", "")
            reordered = self.reorder("p" + out_order).data()
            
            for t in range(self.nv / 2):
                tag = 2*t
                ctrl = tag+1
                output_data[:,:,:,t] = reordered[:,:,:,ctrl] - reordered[:,:,:,tag]
        
        out_order = self.order.replace("p", "").replace("P", "")
        return AslImage(self.ipath + "_diff", data=output_data, 
                        order=out_order, tis=self.tis, ntis=self.ntis, rpts=self.rpts, base=self)

    def mean_across_repeats(self):
        if self.ntc == 2:
            # Have tag-control pairs - need to diff
            diff = self.diff()
        else:
            diff = self
        
        # Reorder so repeats are together
        diff = diff.reorder("tr")
        input_data = diff.data()

        # Create output data - one volume per ti
        output_data = np.zeros(list(self.shape[:3]) + [self.ntis])
        start = 0
        for ti, nrp in enumerate(self.rpts):
            repeat_data = input_data[:,:,:,start:start+nrp]
            output_data[:,:,:,ti] = np.mean(repeat_data, 3)
            start += nrp
        
        return AslImage(self.ipath + "_mean", data=output_data, 
                       order=diff.order, tis=self.tis, ntis=self.ntis, nrpts=1,
                       base=self)

    def summary(self, log=sys.stdout):
        ti_str = "TIs "
        if self.plds: ti_str = "PLDs"
        fsl.Image.summary(self, log)
        log.write("Data shape                    : %s\n" % str(self.shape))
        log.write("Number of %s                : %i\n" % (ti_str, self.ntis))
        log.write("Number of repeats at each TI  : %s\n" % str(self.rpts))
        log.write("Label-Control                 : ")
        if self.ntc == 2:
            if self.tagfirst: log.write("Label-control pairs\n")
            else: log.write("Control-Label pairs\n")
        else:
            log.write("Already differenced\n")

    def derived(self, data, name=None, suffix=None, **kwargs):
        """
        Create a derived ASL image based on this one, but with different data

        This is only possible if the number of volumes match, otherwise we cannot
        use the existing information about TIs, repeats etc. If the number of volumes
        do not match a generic fsl.Image is returned instead
        
        :param data: Numpy data for derived image
        :param name: Name for new image (can be simple name or full filename)
        :param suffix: If name not specified, construct by adding suffix to original image name

        Any further keyword parameters are passed to the Image constructor
        """
        if data.ndim != 4 or data.shape[3] != self.shape[3]:
            return fsl.Image.derived(self, data, name=name, suffix=suffix, **kwargs)
        else:
            if name is None and suffix is None:
                name = self.ipath
            elif name is None:
                name = self.ipath + suffix
            return AslImage(name, data=data, base=self,
                            order=self.order, ntis=self.ntis, tis=self.tis, rpts=self.rpts, **kwargs)

class AslWorkspace(fsl.Workspace):
        
    def smooth(self, img, fwhm, output_name=None):
        if output_name is None:
            output_name = img.iname + "_smooth"
        sigma = round(fwhm/2.355, 2)
        self.log.write("Spatial smoothing with FWHM: %f (sigma=%f)\n" % (fwhm, sigma))
        args = "%s -kernel gauss %f -fmean %s" % (img.iname, sigma, output_name)
        imgs, files, stdout = self.run("fslmaths", args=args, expected=[output_name])
        return imgs[0]

    def preprocess(self, asldata, options, ref=None):
        self.log.write("ASL preprocessing...\n")

        if options.diff: 
            self.log.write("Tag-control subtraction\n")
            asldata = asldata.diff().reorder("rt").save()

        # Keep original AslImage with info about TIs, repeats, etc
        orig = asldata
        if options.mc: 
            self.log.write("Motion correction\n")
            asldata = self.mcflirt(asldata, ref=ref, cost="mutualinfo")

        if options.smooth: 
            asldata = self.smooth(asldata, fwhm=options.fwhm)

        self.log.write("DONE\n\n")
        return orig.derived(data=asldata.data(), name=asldata.ipath)

    def reg(wsp, ref, reg_targets, options, ref_str="asl"):
        """ FIXME not functional yet"""
        self.log.write("Segmentation and co-registration...\n")

        # Brain-extract ref image
        ref_bet = wsp.bet(ref)

        # This is done to avoid the contrast enhanced rim resulting from low intensity ref image
        d = ref_bet.data()
        thr_ref = np.percentile(d[d!=0], 10.0)
        d[d<thr_ref] = 0
        raw_bet = ref_bet.derived(d, save=True)

        for imgs in reg_targets:
            reg = imgs[0]
            reg_bet = wsp.bet(reg, args="-B -f 0.3")
            name = "%s_2%s" % (reg.iname, ref_str)
            name_inv = "%s_2%s" % (ref_str, reg.iname)
            postreg, mat, invmat = wsp.flirt(ref_bet, args="-dof 7", 
                                             output_name=name,
                                             output_mat=name + ".mat",
                                             output_invmat=name_inv + ".mat")
            wsp.write_file(matrix_to_text(invmat), name_inv)
            for coreg in imgs[1:]:
                wsp.apply_xfm(coreg, ref_bet, name_inv, args="-interp nearestneighbour", output_name="%s_fast_seg_2asl" % t1.iname)
    
        self.log.write("DONE\n\n")

def main():
    try:
        p = AslOptionParser(usage="asl_preproc", version="@VERSION@")
        p.add_option("-o", dest="output", help="Output name", default=None)
        options, args = p.parse_args(sys.argv)

        data = AslImage(options.name, **vars(options))
        if options.output is None:
            options.output = data.iname + "_out"
        data.summary()
        
        wsp = AslWorkspace()
        data_preproc = wsp.preprocess(asldata, options)
        data_preproc.save(options.output)

    except Exception as e:
        print("ERROR: " + str(e) + "\n")
        raise

if __name__ == "__main__":
    main()

