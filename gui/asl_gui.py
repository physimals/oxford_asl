#!/usr/bin/env python
import os
import colorsys

import wx
import wx.grid

class RunSequence:
    def __init__(self):
        self.cmds = []

    def add(self, cmd):
        self.cmds.append(cmd)

    def __str__(self):
        s = ""
        for cmd in self.cmds: s += str(cmd) + "\n\n"
        return s

class AslCmd():
    def __init__(self, cmd):
        self.cmd = cmd
    
    def add(self, opt, val=None):
        if val is not None:
            self.cmd += " %s=%s" % (opt, str(val))
        else:
            self.cmd += " %s" % opt

    def __str__(self): return self.cmd

class TabPage(wx.Panel):

    def __init__(self, parent, title, name=None):
        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)
 
        self.sizer = wx.GridBagSizer(vgap=5, hgap=5)
        self.row = 0
        self.title = title
        if name is None:
            self.name = title.lower()
        else:
            self.name = name
            
    def pack(self, label, *widgets, **kwargs):
        col = 0
        border = kwargs.get("border", 10)
        font = self.GetFont()
        if "size" in kwargs:
            font.SetPointSize(kwargs["size"])
        if kwargs.get("bold", False):
            font.SetWeight(wx.BOLD)
        
        if label != "":
            text = wx.StaticText(self, label=label)
            text.SetFont(font)
            self.sizer.Add(text, pos=(self.row, col), border=border, flag=wx.ALIGN_CENTRE_VERTICAL | wx.LEFT)
            col += 1
        else:
            text = None

        for w in widgets:
            span = (1, 1)
            w.label = text
            if hasattr(w, "span"): span = (1, w.span)
            w.SetFont(font)
            w.Enable(col == 0 or kwargs.get("enable", True))
            self.sizer.Add(w, pos=(self.row, col), border=border, flag=wx.ALIGN_CENTRE_VERTICAL | wx.EXPAND | wx.LEFT, span=span)
            col += span[1]
        self.row += 1

    def file_picker(self, label, dir=False, handler=None, optional=False, initial_on=False, pack=True, **kwargs):
        if not handler: handler = self.update
        if dir: 
            picker = wx.DirPickerCtrl(self)
            picker.Bind(wx.EVT_DIRPICKER_CHANGED, handler)
        else: 
            picker = wx.FilePickerCtrl(self)
            picker.Bind(wx.EVT_FILEPICKER_CHANGED, handler)
        picker.span = 2
        if optional:
            cb = wx.CheckBox(self, label=label)
            cb.SetValue(initial_on)
            cb.Bind(wx.EVT_CHECKBOX, handler)
            picker.checkbox = cb
            if pack: self.pack("", cb, picker, enable=initial_on, **kwargs)
        elif pack:
            self.pack(label, picker, **kwargs)

        return picker

    def choice(self, label, choices, initial=0, optional=False, initial_on=False, handler=None, pack=True, **kwargs):
        if not handler: handler = self.update
        ch = wx.Choice(self, choices=choices)
        ch.SetSelection(initial)
        ch.Bind(wx.EVT_CHOICE, handler)
        if optional:
            cb = wx.CheckBox(self, label=label)
            cb.SetValue(initial_on)
            cb.Bind(wx.EVT_CHECKBOX, self.update)
            #cb.Bind(wx.EVT_CHECKBOX, self.enabler(ch))
            ch.checkbox = cb
            if pack: self.pack("", cb, ch, enable=initial_on, **kwargs)
        elif pack:
            self.pack(label, ch, **kwargs)
        return ch

    def number(self, label, handler=None, **kwargs):
        if not handler: handler = self.update
        num = NumberChooser(self, changed_handler=handler, **kwargs)
        num.span = 2
        self.pack(label, num, **kwargs)
        return num

    def integer(self, label, handler=None, pack=True, **kwargs):
        if not handler: handler = self.update
        spin = wx.SpinCtrl(self, **kwargs)
        spin.Bind(wx.EVT_SPINCTRL, handler)
        if pack: self.pack(label, spin)
        return spin

    def checkbox(self, label, initial=False, handler=None, **kwargs):
        cb = wx.CheckBox(self, label=label)
        cb.span=2
        cb.SetValue(initial)
        if handler: cb.Bind(wx.EVT_CHECKBOX, handler)
        else: cb.Bind(wx.EVT_CHECKBOX, self.update)
        self.pack("", cb, **kwargs)
        return cb

    def section(self, label):
        self.pack(label, bold=True)

    def enabler(self, w):
        def enable(evt):
            w.Enable(evt.IsChecked())
        return enable

    def update(self, evt=None):
        if hasattr(self, "run"): self.run.update()

class AslRun(wx.Frame):

    order_opts = {"trp" : "--ibf=tis --iaf=diff", 
                  "trp,tc" : "--ibf=tis --iaf=tcb --diff", 
                  "trp,ct" : "--ibf=tis --iaf=ctb --diff",
                  "rtp" : "--ibf=rpt --iaf=diff",
                  "rtp,tc" : "--rpt --iaf=tcb --diff",
                  "rtp,ct" : "--ibf=rpt --iaf=ctb --diff",
                  "ptr,tc" : "--ibf=tis --iaf=tc --diff",
                  "ptr,ct" : "--ibf=tis --iaf=ct --diff",
                  "prt,tc" : "--ibf=rpt --iaf=tc --diff",
                  "prt,ct" : "--ibf=rpt --iaf=ct --diff"}

    def __init__(self, parent, run_btn, run_label):
        wx.Frame.__init__(self, parent, title="Run", size=(600, 400), style=wx.DEFAULT_FRAME_STYLE)

        self.run_btn = run_btn
        self.run_label = run_label
        self.run_btn.Bind(wx.EVT_BUTTON, self.dorun)

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.output_text = wx.TextCtrl(self, style=wx.TE_READONLY | wx.TE_MULTILINE)
        self.sizer.Add(self.output_text, 1, flag=wx.EXPAND)
            
        self.SetSizer(self.sizer)
        self.Bind(wx.EVT_CLOSE, self.close)

    def close(self, evt):
        self.Hide()

    def dorun(self, evt):
        if self.run_seq: 
            self.output_text.Clear()
            self.output_text.AppendText(str(self.run_seq))
            self.Show()

    def update(self, evt=None):
        self.run_seq = None
        try:
            self.run_seq = self.get_run_sequence()
            self.run_label.SetForegroundColour(wx.Colour(0, 255, 0))
            self.run_label.SetLabel("Ready to Go")
            self.run_btn.Enable(True)
        except Exception, e:
            self.run_btn.Enable(False)
            self.run_label.SetForegroundColour(wx.Colour(255, 0, 0))
            self.run_label.SetLabel(str(e))

    def check_exists(self, label, file):
        if not os.path.exists(file):
            raise RuntimeError("%s - no such file or directory" % label)

    def get_run_sequence(self):
        run = RunSequence()

        self.check_exists("Input data", self.input.data())

        # Create output dirs
        outdir = self.analysis.outdir()
        if outdir == "": 
            raise RuntimeError("No output dir")
        run.add("mkdir %s" % outdir)
        run.add("mkdir %s/native_space" % outdir)

        # Check data order is supported
        order, tagfirst = self.input.data_order()
        if self.input.tc_pairs(): 
            if tagfirst: order += ",tc"
            else: order += ",ct"
        if order not in self.order_opts:
            #print(order)
            raise RuntimeError("This data ordering is not supported by ASL_FILE")
        
        # OXFORD_ASL
        cmd = AslCmd("oxford_asl")
        cmd.add("-i %s" % self.input.data())
        cmd.add("-o %s" % outdir)
        cmd.add(self.order_opts[order])
        cmd.add("--tis %s" % ",".join(["%.2f" % v for v in self.input.tis()]))
        cmd.add("--bolus %s" % ",".join(["%.2f" % v for v in self.input.bolus_dur()]))
        if self.analysis.wp(): 
            cmd.add("--wp")
        else: 
            cmd.add("--t1 %.2f" % self.analysis.t1())
            cmd.add("--bat %.2f" % self.analysis.bat())
        cmd.add("--t1b %.2f" % self.analysis.t1b())
        cmd.add("--alpha %.2f" % self.analysis.ie())
        if self.analysis.spatial(): cmd.add("--spatial")
        if self.analysis.infer_t1(): cmd.add("--infert1")
        if self.analysis.fixbolus(): cmd.add("--fixbolus")
        if self.analysis.pv(): cmd.add("--pvcorr")
        if self.analysis.mc(): cmd.add("--mc")
        if not self.analysis.macro(): cmd.add("--artoff")
        if self.analysis.mask() is not None:
            self.check_exists("Analysis mask", self.analysis.mask())
            cmd.add("-m %s" % self.analysis.mask())
        if self.input.labelling() == 1: 
            cmd.add("--casl")
        if self.analysis.transform():
            if self.analysis.transform_type() == 0:
                self.check_exists("Transformation matrix", self.analysis.transform_file())
                cmd.add("--asl2struc %s" % self.analysis.transform_file())
            elif self.analysis.transform_type() == 1:
                self.check_exists("Warp image", self.analysis.transform_file())
                cmd.add("--regfrom %s" % self.analysis.transform_file())
            else:
                pass # --fslanat already set when option 2 chosen

        # Structural image - may require Bet to be run
        fsl_anat_dir = self.input.fsl_anat()
        struc_image = self.input.struc_image()
        if fsl_anat_dir is not None:
            self.check_exists("FSL_ANAT", fsl_anat_dir)
            cmd.add("--fslanat=%s" % fsl_anat_dir)
        elif struc_image is not None:
            self.check_exists("Structural image", struc_image)
            cp = AslCmd("imcp")
            cp.add(struc_image)
            cp.add("%s/structural_head" % outdir)
            run.add(cp)
            if self.input.struc_image_bet() == 1:
                bet = AslCmd("bet")
                bet.add(struc_image)
                bet.add("%s/structural_brain" % outdir)
                run.add(bet)
            else:
                struc_image_brain = self.input.struc_image_brain()
                self.check_exists("Structural brain image", struc_image_brain)
                cp = AslCmd("imcp")
                cp.add(struc_image_brain)
                cp.add("%s/structural_brain" % outdir)
                run.add(cp)
            cmd.add("--s %s/structural_head" % outdir)
            cmd.add("--sbrain %s/structural_brain" % outdir)
        else:
            # No structural image
            pass
        
        run.add(cmd)
    
        # ASL_CALIB
        if self.calibration.calib():
            calib = AslCmd("asl_calib")
            calib.add("-i %s/native_space/perfusion" % outdir)
            calib.add("-o %s/calibration" % outdir)
            calib.add("-c %s" % self.calib.calib_image())
            if self.calib.m0_type() == 0:
                calib.add("--mode longtr")
                calib.add("--tr %.2f" % self.calib.seq_tr())
            else:
                calib.add("--mode satrevoc")
                calib.add("--tis %s" % ",".join([str(v) for v in self.input.tis()]))
                # FIXME change -c option in sat recov mode?

            # FIXME structural image required?
            #calib.add("-s %s" % self.calib.calib_image())
            calib.add("--te %.2f" % self.calibration.seq_te())
            if self.calibration.calib_mode() == 0:
                calib.add("--tissref %s" % self.calib.ref_tissue_type_name())
                calib.add("--t1r %.2f" % self.calibration.ref_t1())
                calib.add("--t2r %.2f" % self.calibration.ref_t2())
                calib.add("--t2b %.2f" % self.calibration.blood_t2())
                if self.calibration.ref_tissue_mask() is not None:
                    self.check_exists("Calibration reference tissue mask", self.calibration.ref_tissue_mask())
                    calib.add("-m %s" % self.calibration.ref_tissue_mask())
                else:
                    # use structural_brain? What about FSL_ANAT?
                    calib.add("-s %s/structural_brain" % self.outdir)
                    calib.add("-t %s/native_space/asl2struct.mat" % self.outdir)
            if self.calibration.coil_image() is not None:
                self.check_exists("Coil sensitivity reference image", self.calibration.coil_image())
                calib.add("--cref %s" % self.calibration.coil_image())

            run.add("FIXME need to add more calibration commands here")

        return run

class AslCalibration(TabPage):

    def __init__(self, parent):
        TabPage.__init__(self, parent, "Calibration")

        self.wp = False

        self.calib_cb = self.checkbox("Enable Calibration", bold=True)

        self.calib_image_picker = self.file_picker("Calibration Image")
        self.m0_type_ch = self.choice("M0 Type", choices=["Proton Density (long TR)", "Saturation Recovery"])

        self.seq_tr_num = self.number("Sequence TR (s)", min=0,max=10,initial=6)
        self.calib_gain_num = self.number("Calibration Gain", min=0,max=5,initial=1)
        self.calib_mode_ch = self.choice("Calibration mode", choices=["Reference Region", "Voxelwise"])

        self.section("Reference tissue")

        self.ref_tissue_type_ch = self.choice("Type", choices=["CSF", "WM", "GM", "None"], handler=self.ref_tissue_type_changed)
        self.ref_tissue_mask_picker = self.file_picker("Mask", optional=True)
        self.seq_te_num = self.number("Sequence TE (ms)", min=0,max=30,initial=0)
        self.blood_t2_num = self.number("Blood T2 (s)", min=0,max=5,initial=3)
        self.coil_image_picker = self.file_picker("Coil Sensitivity Image", optional=True)
        self.ref_t1_num = self.number("Reference T1 (s)", min=0,max=5,initial=1.3)
        self.ref_t2_num = self.number("Reference T2 (s)", min=0,max=5,initial=1)

        self.sizer.AddGrowableCol(2, 1)
        #self.sizer.AddGrowableRow(self.row, 1)
        self.SetSizer(self.sizer)

    def calib(self): return self.calib_cb.IsChecked()
    def m0_type(self): return self.m0_type_ch.GetSelection()
    def seq_tr(self): return self.seq_tr_num.GetValue()
    def seq_te(self): return self.seq_te_num.GetValue()
    def calib_image(self): return self.calib_image_picker.GetPath()
    def calib_gain(self): return self.calib_gain_num.GetValue()
    def calib_mode(self): return self.calib_mode_ch.GetSelection()
    def ref_tissue_type(self): return self.ref_tissue_type_ch.GetSelection()
    def ref_tissue_type_name(self): return self.ref_tissue_type_ch.GetString(self.ref_tissue_type())
    def ref_tissue_mask(self): 
        if self.ref_tissue_mask_picker.checkbox.IsChecked():
            return self.ref_tissue_mask_picker.GetPath()
        else:
            return None
    def ref_t1(self): return self.ref_t1_num.getValue()
    def ref_t2(self): return self.ref_t2_num.getValue()
    def blood_t2(self): return self.blood_t2_num.getValue()
    def coil_image(self): 
        if self.coil_image_picker.checkbox.IsChecked: return self.coil_image_picker.GetPath()
        else: return None

    def ref_tissue_type_changed(self, event):
        if self.ref_tissue_type() == 0: # CSF
            self.ref_t1_num.SetValue(4.3)
            self.ref_t2_num.SetValue(750.0/400)
        elif self.ref_tissue_type() == 1: # WM
            self.ref_t1_num.SetValue(1.0)
            self.ref_t2_num.SetValue(50.0/50)
        elif self.ref_tissue_type() == 2: # GM
            self.ref_t1_num.SetValue(1.3)
            self.ref_t2_num.SetValue(100.0/60)
        self.update()

    def wp_changed(self, wp):
        self.wp = wp
        if wp: self.calib_mode_ch.SetSelection(1)
        self.update()

    def update(self, event=None):
        enable = self.calib()
        self.m0_type_ch.Enable(enable)
        self.seq_tr_num.Enable(enable and self.m0_type() == 0)
        self.calib_image_picker.Enable(enable)
        self.calib_gain_num.Enable(enable)
        self.coil_image_picker.checkbox.Enable(enable)
        self.calib_mode_ch.Enable(enable and not self.wp)
        self.ref_tissue_type_ch.Enable(enable and self.calib_mode() == 0)
        self.ref_tissue_mask_picker.checkbox.Enable(enable and self.calib_mode() == 0)
        self.ref_tissue_mask_picker.Enable(enable and self.ref_tissue_mask_picker.checkbox.IsChecked() and self.calib_mode() == 0)
        self.coil_image_picker.checkbox.Enable(enable and self.calib_mode() == 0)
        self.coil_image_picker.Enable(enable and self.calib_mode() == 0 and self.coil_image_picker.checkbox.IsChecked())
        self.seq_te_num.Enable(enable and self.calib_mode() == 0)
        self.blood_t2_num.Enable(enable and self.calib_mode() == 0)
        self.ref_t1_num.Enable(enable and self.calib_mode() == 0)
        self.ref_t2_num.Enable(enable and self.calib_mode() == 0)
        TabPage.update(self)

class AslAnalysis(TabPage):

    def __init__(self, parent):
        TabPage.__init__(self, parent, "Analysis")

        self.transform_choices = ["Matrix", "Warp image", "Use FSL_ANAT output"]

        self.section("Registration")

        self.transform_cb = wx.CheckBox(self, label="Transform to standard space")
        self.transform_cb.Bind(wx.EVT_CHECKBOX, self.update)
        self.transform_ch = wx.Choice(self, choices=self.transform_choices)
        self.transform_ch.SetSelection(2)
        self.transform_ch.Bind(wx.EVT_CHOICE, self.update)
        self.transform_picker = wx.FilePickerCtrl(self)
        self.pack("", self.transform_cb, self.transform_ch, self.transform_picker, enable=False)

        self.section("Basic analysis options")

        self.outdir_picker = self.file_picker("Output Directory", dir=True)
        self.mask_picker = self.file_picker("Brain Mask", optional=True)
        self.wp_cb = self.checkbox("Analysis which conforms to 'White Paper' (Alsop et al 2014)", handler=self.wp_changed)

        self.section("Initial parameter values")

        self.bat_num = self.number("Bolus arrival time (s)", min=0,max=2.5,initial=0.7)
        self.t1_num = self.number("T1 (s)", min=0,max=3,initial=1.3)
        self.t1b_num = self.number("T1b (s)", min=0,max=3,initial=1.65)
        self.ie_num = self.number("Inversion Efficiency", min=0,max=1,initial=0.98)
        
        self.section("Analysis Options")

        self.spatial_cb = self.checkbox("Adaptive spatial regularization on perfusion", initial=True)
        self.infer_t1_cb = self.checkbox("Incorporate T1 value uncertainty")
        self.macro_cb = self.checkbox("Include macro vascular component")
        self.fixbolus_cb = self.checkbox("Fix bolus duration", initial=True)

        self.pv_cb = self.checkbox("Partial Volume Correction")
        self.mc_cb = self.checkbox("Motion Correction (MCFLIRT)")

        self.sizer.AddGrowableCol(1, 1)
        #sizer.AddGrowableRow(5, 1)
        self.SetSizer(self.sizer)

    def transform(self): return self.transform_cb.IsChecked()
    def transform_type(self): return self.transform_ch.GetSelection()
    def transform_file(self): return self.transform_picker.GetPath()
    def outdir(self): return self.outdir_picker.GetPath()
    def mask(self): 
        if self.mask_picker.checkbox.IsChecked(): return self.mask_picker.GetPath()
        else: return None
    def wp(self): return self.wp_cb.IsChecked()
    def bat(self): return self.bat_num.GetValue()
    def t1(self): return self.t1_num.GetValue()
    def t1b(self): return self.t1b_num.GetValue()
    def ie(self): return self.ie_num.GetValue()
    def spatial(self): return self.spatial_cb.IsChecked()
    def infer_t1(self): return self.infer_t1_cb.IsChecked()
    def macro(self): return self.macro_cb.IsChecked()
    def fixbolus(self): return self.fixbolus_cb.IsChecked()
    def pv(self): return self.pv_cb.IsChecked()
    def mc(self): return self.mc_cb.IsChecked()

    def update(self, event=None):
        self.transform_ch.Enable(self.transform())
        self.transform_picker.Enable(self.transform() and self.transform_type() != 2)
        self.mask_picker.Enable(self.mask_picker.checkbox.IsChecked())
        TabPage.update(self)

    def wp_changed(self, event):
        if self.wp():
            self.t1_num.SetValue(1.65)
        self.calibration.wp_changed(self.wp())
        self.t1_num.Enable(not self.wp())
        self.update()

    def labelling_changed(self, pasl):
        if pasl:
            self.bat_num.SetValue(0.7)
            self.ie_num.SetValue(0.98)
        else:
            self.bat_num.SetValue(1.3)
            self.ie_num.SetValue(0.85)
        
    def fsl_anat_changed(self, enabled):
        """ If FSL_ANAT is selected, use it by default, otherwise do not allow it """
        sel = self.transform_ch.GetSelection()
        if enabled: 
            choices = self.transform_choices
            sel = 2
        else: 
            choices = self.transform_choices[:2]
            if sel == 2: sel = 0
        self.transform_ch.Enable(False)
        self.transform_ch.Clear()
        self.transform_ch.AppendItems(choices)
        self.transform_ch.SetSelection(sel)
        self.transform_ch.Enable(self.transform())

class AslInputOptions(TabPage):
    
    def __init__(self, parent):
        TabPage.__init__(self, parent, "Input Data", "input")
 
        self.groups = ["TIs", "Repeats", "Tag/Control pairs"]
        self.abbrevs = ["t", "r", "p"]

        self.section("Data contents")

        self.data_picker = self.file_picker("Input Image")
        self.ntis_int = self.integer("Number of TIs", min=1,max=100,initial=1)
        self.nrepeats_int = self.integer("Number of repeats", min=1,max=100,initial=1)

        self.section("Data order")

        self.choice1 = wx.Choice(self, choices=self.groups)
        self.choice1.SetSelection(2)
        self.choice1.Bind(wx.EVT_CHOICE, self.update)
        self.choice2 = wx.Choice(self, choices=self.groups)
        self.choice2.SetSelection(0)
        self.choice2.Bind(wx.EVT_CHOICE, self.update)
        self.pack("Grouping order", self.choice1, self.choice2)
        self.tc_ch = self.choice("Tag control pairs", choices=["Tag then control", "Control then tag"], optional=True, initial_on=True)
        
        self.section("Acquisition parameters")

        self.labelling_ch = self.choice("Labelling", choices=["pASL", "cASL/pcASL"], initial=1, handler=self.labelling_changed)

        self.bolus_dur_ch = wx.Choice(self, choices=["Constant", "Variable"])
        self.bolus_dur_ch.SetSelection(0)
        self.bolus_dur_ch.Bind(wx.EVT_CHOICE, self.update)
        self.bolus_dur_num = NumberChooser(self, min=0, max=2.5, step=0.1, initial=0.7)
        self.bolus_dur_num.span = 2
        self.pack("Bolus duration (s)", self.bolus_dur_ch, self.bolus_dur_num)
        
        self.bolus_dur_list = NumberList(self, self.ntis())
        self.bolus_dur_list.span = 3
        self.bolus_dur_list.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self.update)
        self.pack("Bolus durations (s)", self.bolus_dur_list, enable=False)

        self.ti_list = NumberList(self, self.ntis())
        self.ti_list.span=3
        self.ti_list.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self.update)
        self.pack("TIs (s)", self.ti_list)
        
        self.readout_ch = wx.Choice(self, choices=["3D (eg GRASE)", "2D multi-slice (eg EPI)"])
        self.readout_ch.SetSelection(0)
        self.readout_ch.Bind(wx.EVT_CHOICE, self.update)
        self.time_per_slice_num = NumberChooser(self, label="Time per slice (ms)", min=0, max=50, step=1, initial=10)
        self.time_per_slice_num.span=2
        self.pack("Readout", self.readout_ch, self.time_per_slice_num)
        self.time_per_slice_num.Enable(False)
        
        self.multiband_cb = wx.CheckBox(self, label="Multi-band")
        self.multiband_cb.Bind(wx.EVT_CHECKBOX, self.update)
        self.slices_per_band_spin = wx.SpinCtrl(self, min=1, max=100, initial=5)
        self.slices_per_band_label = wx.StaticText(self, label="slices per band")
        self.pack("", self.multiband_cb, self.slices_per_band_spin, self.slices_per_band_label, enable=False)
        self.multiband_cb.Enable(False)

        self.section("Structure")

        self.fsl_anat_picker = self.file_picker("Use FSL_ANAT output", dir=True, optional=True, initial_on=True, handler=self.fsl_anat_changed)
        self.struc_image_picker = self.file_picker("Use Structural Image", optional=True)
        self.struc_image_ch = wx.Choice(self, choices=["Have brain image", "Run BET to extract brain"])
        self.struc_image_ch.Bind(wx.EVT_CHOICE, self.update)
        self.struc_image_ch.SetSelection(0)
        self.struc_image_brain_picker = wx.FilePickerCtrl(self)
        self.struc_image_brain_picker.span = 2
        self.struc_image_brain_picker.Bind(wx.EVT_FILEPICKER_CHANGED, self.update)
        self.pack("Brain extraction", self.struc_image_ch, self.struc_image_brain_picker)

        self.section("Data Order Preview")

        self.preview = AslDataPreview(self, self.ntis(), self.nrepeats(), self.tc_pairs(), "trp", True)
        self.sizer.Add(self.preview, pos=(self.row, 0), span=(1, 5), flag=wx.EXPAND | wx.ALL)

        self.sizer.AddGrowableCol(2, 1)
        self.sizer.AddGrowableRow(self.row, 1)
        self.SetSizer(self.sizer)

    def data(self): return self.data_picker.GetPath()
    def ntis(self): return self.ntis_int.GetValue()
    def nrepeats(self): return self.nrepeats_int.GetValue()
    def data_order(self): return self.preview.order, self.preview.tagfirst
    def tc_pairs(self): return self.tc_ch.checkbox.IsChecked()
    def labelling(self): return self.labelling_ch.GetSelection()
    def bolus_dur_type(self): return self.bolus_dur_ch.GetSelection()
    def bolus_dur(self): 
        if self.bolus_dur_type() == 0: return [self.bolus_dur_num.GetValue(), ]
        else: return self.bolus_dur_list.GetValues()
    def tis(self): return self.ti_list.GetValues()
    def readout(self): return self.readout_ch.GetSelection()
    def time_per_slice(self): return self.time_per_slice_num.GetValue()
    def multiband(self): return self.multiband_cb.IsChecked()
    def slices_per_band(self): return self.slices_per_band_spin.GetValue()
    def fsl_anat(self): 
        if self.fsl_anat_picker.checkbox.IsChecked(): return self.fsl_anat_picker.GetPath()
        else: return None
    def struc_image(self): 
        if self.struc_image_picker.checkbox.IsChecked() and not self.fsl_anat():
            return self.struc_image_picker.GetPath()
        else: return None
    def struc_image_brain(self): 
        if self.struc_image_picker.checkbox.IsChecked() and not self.fsl_anat() and self.struc_image_bet() == 0:
            return self.struc_image_brain_picker.GetPath()
        else: return None
    def struc_image_bet(self): return self.struc_image_ch.GetSelection()

    def update(self, event=None):
        self.ti_list.set_size(self.ntis())
        self.bolus_dur_list.set_size(self.ntis())

        self.time_per_slice_num.Enable(self.readout() != 0)
        self.multiband_cb.Enable(self.readout() != 0)
        self.slices_per_band_spin.Enable(self.multiband() and self.readout() != 0)
        self.slices_per_band_label.Enable(self.multiband() and self.readout() != 0)

        use_fsl_anat = self.fsl_anat_picker.checkbox.IsChecked()
        use_struc_image = self.struc_image_picker.checkbox.IsChecked()
        self.fsl_anat_picker.Enable(use_fsl_anat)
        self.struc_image_picker.checkbox.Enable(not use_fsl_anat)
        self.struc_image_ch.Enable(use_struc_image and not use_fsl_anat)
        self.struc_image_picker.Enable(use_struc_image and not use_fsl_anat)
        self.struc_image_brain_picker.Enable(use_struc_image and not use_fsl_anat and self.struc_image_bet() == 0)
        self.struc_image_brain_picker.label.Enable(use_struc_image and not use_fsl_anat and self.struc_image_bet() == 0)
        
        self.bolus_dur_num.Enable(self.bolus_dur_type() == 0)
        self.bolus_dur_list.Enable(self.bolus_dur_type() == 1)

        self.tc_ch.Enable(self.tc_pairs())
        self.update_groups()

        self.preview.n_tis = self.ntis()
        self.preview.n_repeats = self.nrepeats()
        self.preview.tc_pairs = self.tc_pairs()
        self.preview.tagfirst = self.tc_ch.GetSelection() == 0
        self.preview.Refresh()
        TabPage.update(self)

    def labelling_changed(self, event):
        if event.GetInt() == 0:
            self.bolus_dur_num.SetValue(0.7)
            self.ntis_int.label.SetLabel("Number of TIs")
            self.ti_list.label.SetLabel("TIs")
        else:
            self.bolus_dur_num.SetValue(1.8)
            self.ntis_int.label.SetLabel("Number of PLDs")
            self.ti_list.label.SetLabel("PLDs")
        self.analysis.labelling_changed(event.GetInt() == 0)
        self.update()

    def fsl_anat_changed(self, event):
        self.analysis.fsl_anat_changed(self.fsl_anat_picker.checkbox.IsChecked())
        self.update()

    def update_groups(self, group1=True, group2=True):
        g2 = self.choice2.GetSelection()
        g1 = self.choice1.GetSelection()
        if not self.tc_pairs():
            if g1 == 2: g1 = 0
            if g2 == 1: g2 = 0
            choices = 2
        else:
            choices = 3
        group1_items = []
        group2_items = []
        for idx, item in enumerate(self.groups[:choices]):
            group1_items.append(item)
            if idx != g1: group2_items.append(item)

        self.update_group_choice(self.choice1, group1_items, g1)
        self.update_group_choice(self.choice2, group2_items, g2)
        
        if g2 >= g1: g2 += 1
        order = self.abbrevs[g1]
        order += self.abbrevs[g2]
        order += self.abbrevs[3-g1-g2]
        #print(order)
        self.preview.order = order

    def update_group_choice(self, w, items, sel):
        w.Enable(False)
        w.Clear()
        w.AppendItems(items)
        w.SetSelection(sel)
        w.Enable(True)

class NumberChooser(wx.Panel):
    def __init__(self, parent, label=None, min=0, max=1, initial=0.5, step=0.1, digits=2, changed_handler=None):
        super(NumberChooser, self).__init__(parent)
        self.min = min
        self.max = max
        self.handler = changed_handler
        self.hbox=wx.BoxSizer(wx.HORIZONTAL)
        if label is not None:
            self.label = wx.StaticText(self, label=label)
            self.hbox.Add(self.label, proportion=0, flag=wx.ALIGN_CENTRE_VERTICAL)
        self.spin = wx.SpinCtrlDouble(self, min=min,max=max,initial=initial)
        self.spin.SetDigits(digits)
        self.spin.SetIncrement(step)
        self.spin.Bind(wx.EVT_SLIDER, self.spin_changed)
        self.slider = wx.Slider(self, value=initial, minValue=0, maxValue=100)
        self.slider.SetValue(100*(initial-self.min)/(self.max-self.min))
        self.slider.Bind(wx.EVT_SLIDER, self.slider_changed)
        self.hbox.Add(self.slider, proportion=1, flag=wx.EXPAND | wx.ALIGN_CENTRE_VERTICAL)
        self.hbox.Add(self.spin, proportion=0, flag=wx.EXPAND | wx.ALIGN_CENTRE_VERTICAL)
        self.SetSizer(self.hbox)

    def GetValue(self):
        return self.spin.GetValue()
        
    def SetValue(self, val):
        self.spin.SetValue(val)
        self.slider.SetValue(100*(val-self.min)/(self.max-self.min))
        
    def slider_changed(self, event):
        v = event.GetInt()
        val = self.min + (self.max-self.min)*float(v)/100
        self.spin.SetValue(val)
        if self.handler: self.handler(event)

    def spin_changed(self, event):
        val = event.GetValue()
        self.slider.SetValue(100*(val-self.min)/(self.max-self.min))
        if self.handler: self.handler(event)

class NumberList(wx.grid.Grid):
    def __init__(self, parent, n, default=1.8):
        super(NumberList, self).__init__(parent, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, 0 )
        self.n=0
        self.default = default
        self.CreateGrid(1, 0)
        self.SetRowLabelSize(0)
        self.SetColLabelSize(0)
        self.set_size(n)
        self.Bind(wx.EVT_SIZE, self.on_size)
        #for i in range(self.n):
        #    self.SetColSize(i, 20)

    def GetValues(self):
        try:
            return [float(self.GetCellValue(0, c)) for c in range(self.n)]
        except ValueError, e:
            raise RuntimeError("Non-numeric values in number list")
            
    def set_size(self, n):
        if n > self.n:
            self.AppendCols(n - self.n)
            for c in range(self.n, n): self.SetCellValue(0, c, str(self.default))
        elif n < self.n:
            self.DeleteCols(n, self.n-n)
        self.n = n
        self.resize_cols()

    def resize_cols(self):
        w, h = self.GetClientSize()
        cw = w / self.n
        for i in range(self.n):
            self.SetColSize(i, cw)

    def on_size(self, event):
        event.Skip()
        self.resize_cols()

class AslDataPreview(wx.Panel):
    def __init__(self, parent, n_tis, n_repeats, tc_pairs, order, tagfirst):
        super(AslDataPreview, self).__init__(parent)
        self.SetBackgroundStyle(wx.BG_STYLE_CUSTOM)
        self.Bind(wx.EVT_SIZE, self.on_size)
        self.Bind(wx.EVT_PAINT, self.on_paint)
        self.n_tis = n_tis
        self.n_repeats = n_repeats
        self.tc_pairs = tc_pairs
        self.tagfirst = tagfirst
        self.order = order

    def on_size(self, event):
        event.Skip()
        self.Refresh()

    def get_col(self, pos, ti):
        if ti: h = 170.0/255
        else: 
            h = 90.0/255
        s, v = 0.5, 0.95 - pos/2
        #print(h, s, v)
        r,g,b = colorsys.hsv_to_rgb(h, s, v)
        #print(r,g,b)
        return wx.Colour(int(r*255), int(g*255), int(b*255))

    def on_paint(self, event):
        w, h = self.GetClientSize()
        N = self.n_tis * self.n_repeats
        if self.tc_pairs: N *= 2
        dc = wx.AutoBufferedPaintDC(self)
        dc.Clear()

        leg_width = (w-200)/4
        leg_start = 100

        #b = wx.Brush(self.get_col(0.5, True), wx.SOLID)
        #dc.SetBrush(b)
        dc.SetBrush(wx.TRANSPARENT_BRUSH)
        rect = wx.Rect(leg_start, 20, leg_width/4, 20)
        dc.GradientFillLinear(rect, self.get_col(0, True), self.get_col(1.0, True), wx.EAST)
        dc.DrawRectangleRect(rect)
        dc.DrawText("Tis", leg_start+leg_width/3, 20)

        #b = wx.Brush(self.get_col(0.5, False), wx.SOLID)
        #dc.SetBrush(b)
        rect = wx.Rect(leg_start+leg_width, 20, leg_width/4, 20)
        dc.GradientFillLinear(rect, self.get_col(0, False), self.get_col(1.0, False), wx.EAST)
        dc.DrawRectangleRect(rect)
        dc.DrawText("Repeats", leg_start+4*leg_width/3, 20)

        if self.tc_pairs:
            dc.SetBrush(wx.TRANSPARENT_BRUSH)
            dc.DrawRectangle(leg_start+leg_width*2, 20, leg_width/4, 20)
            dc.DrawText("Tag", leg_start+7*leg_width/3, 20)

            b = wx.Brush('black', wx.BDIAGONAL_HATCH)
            dc.SetBrush(b)
            dc.DrawRectangle(leg_start+leg_width*3, 20, leg_width/4, 20)
            dc.DrawText("Control", leg_start+10*leg_width/3, 20)

        dc.DrawRectangle(50, 50, w-100, h-100)
        dc.DrawRectangle(50, 50, w-100, h-100)
        dc.DrawText("0", 50, h-50)
        dc.DrawText(str(N), w-50, h-50)

        seq = [1,]
        for t in self.order[::-1]:
            temp = seq
            seq = []
            for i in temp:
                if t == "t":
                    seq += [i,] * self.n_tis
                elif t == "p":
                    if not self.tc_pairs:
                        seq.append(i)
                    elif self.tagfirst:
                        seq.append(i)
                        seq.append(i+1)
                    else:
                        seq.append(i+1)
                        seq.append(i)
                elif t == "r":
                    seq.append(i)
                    seq += [i+2,] * (self.n_repeats - 1)
        
        tistart = -1
        ti_sep = 1
        for idx, s in enumerate(seq):
            if s == 1 and tistart < 0: 
                tistart = idx
            elif s == 1:
                ti_sep = idx - tistart
                break

        bwidth = float(w - 100) / N
        x = 50
        pos = 0.0
        ti = 0
        d = 1.0/self.n_tis
        for idx, s in enumerate(seq):
            dc.SetPen(wx.TRANSPARENT_PEN)
            b = wx.Brush(self.get_col(pos, s in (1, 2)), wx.SOLID)
            dc.SetBrush(b)
            dc.DrawRectangle(int(x), 50, int(bwidth+1), h-100)

            if s in (2, 4):
                b = wx.Brush('black', wx.BDIAGONAL_HATCH)
                dc.SetBrush(b)
                dc.DrawRectangle(int(x), 50, int(bwidth+1), h-100)

            dc.SetPen(wx.Pen('black'))
            dc.DrawLine(int(x), 50, int(x), h-50)
            x += bwidth
            if (idx+1) % ti_sep == 0: 
                pos += d
                ti += 1
                if ti == self.n_tis: 
                    pos = 0
                    ti = 0

class AslGui(wx.Frame):

    def __init__(self):
        wx.Frame.__init__(self, None, title="Basil", size=(800, 850), style=wx.DEFAULT_FRAME_STYLE ^ wx.RESIZE_BORDER)
        panel = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)
        
        banner = wx.StaticBitmap(panel, -1, wx.Bitmap("banner.png", wx.BITMAP_TYPE_ANY))
        sizer.Add(banner)

        notebook = wx.Notebook(panel, id=wx.ID_ANY, style=wx.BK_DEFAULT)
        
        sizer.Add(notebook, 1, wx.ALL|wx.EXPAND, 5)
                
        self.run_panel = wx.Panel(panel)
        runsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.run_label = wx.StaticText(self.run_panel, label="Unchecked")
        self.run_label.SetFont(wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        runsizer.Add(self.run_label, 1, wx.EXPAND)
        self.run_btn = wx.Button(self.run_panel, label="Run")
        runsizer.Add(self.run_btn, 0, wx.ALIGN_CENTER_VERTICAL)
        self.run_panel.SetSizer(runsizer)
        sizer.Add(self.run_panel, 0, wx.EXPAND)
        self.run_panel.Layout()

        self.run = AslRun(self, self.run_btn, self.run_label)
        tabs = [AslInputOptions(notebook),
                AslAnalysis(notebook),
                AslCalibration(notebook)]

        for tab in tabs:
            notebook.AddPage(tab, tab.title)
            setattr(tab, "run", self.run)
            setattr(self.run, tab.name, tab)
            for tab2 in tabs:
                if tab != tab2: setattr(tab, tab2.name, tab2)
            tab.update()

        panel.SetSizer(sizer)
        self.Layout()

if __name__ == '__main__':
    app = wx.App(redirect=False)
    top = AslGui()
    top.Show()
    app.MainLoop()

