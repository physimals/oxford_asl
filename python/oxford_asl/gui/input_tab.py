import sys, os

import wx
import wx.grid

import nibabel as nib

from .widgets import TabPage, NumberChooser, NumberList

class AslInputOptions(TabPage):
    """
    Tab page containing input data options
    """

    def __init__(self, parent, idx, n):
        TabPage.__init__(self, parent, "Input Data", idx, n, name="input")
 
        self.groups = ["PLDs", "Repeats", "Label/Control pairs"]
        self.abbrevs = ["t", "r", "p"]
        self.nvols = -1

        self.section("Data contents")

        self.data_picker = self.file_picker("Input Image", handler=self.set_default_dir)
        self.ntis_int = self.integer("Number of PLDs", min=1,max=100,initial=1)
        self.repeats_choice = self.choice("Repeats", choices=["Fixed", "Variable"], handler=self.update)

        self.section("Data order")

        self.choice1 = wx.Choice(self, choices=self.groups)
        self.choice1.SetSelection(2)
        self.choice1.Bind(wx.EVT_CHOICE, self.update)
        self.choice2 = wx.Choice(self, choices=self.groups)
        self.choice2.SetSelection(0)
        self.choice2.Bind(wx.EVT_CHOICE, self.update)
        self.pack("Grouping order", self.choice1, self.choice2)
        self.tc_ch = self.choice("Label/Control pairs", choices=["Label then control", "Control then label"], optional=True, initial_on=True)
        
        self.section("Acquisition parameters")

        self.labelling_ch = self.choice("Labelling", choices=["pASL", "cASL/pcASL"], initial=1, handler=self.labelling_changed)

        self.bolus_dur_ch = wx.Choice(self, choices=["Constant", "Variable"])
        self.bolus_dur_ch.SetSelection(0)
        self.bolus_dur_ch.Bind(wx.EVT_CHOICE, self.update)
        self.bolus_dur_num = NumberChooser(self, min=0, max=2.5, step=0.1, initial=1.8)
        self.bolus_dur_num.span = 2
        self.bolus_dur_num.spin.Bind(wx.EVT_SPINCTRL, self.bolus_dur_changed)
        self.bolus_dur_num.slider.Bind(wx.EVT_SLIDER, self.bolus_dur_changed)
        self.pack("Bolus duration (s)", self.bolus_dur_ch, self.bolus_dur_num)
        
        self.bolus_dur_list = NumberList(self, self.ntis())
        self.bolus_dur_list.span = 3
        self.bolus_dur_list.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self.update)
        self.pack("Bolus durations (s)", self.bolus_dur_list, enable=False)

        self.ti_list = NumberList(self, self.ntis())
        self.ti_list.span=3
        self.ti_list.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self.update)
        self.pack("PLDs (s)", self.ti_list)
        
        self.nrepeats_list = NumberList(self, self.ntis(), default=1)
        self.nrepeats_list.span = 3
        self.nrepeats_list.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self.update)
        self.pack("Repeats", self.nrepeats_list, enable=False)

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

        self.sizer.AddGrowableCol(2, 1)
        self.SetSizer(self.sizer)
        self.next_prev()

    def data(self): return self.data_picker.GetPath()
    def ntis(self): return self.ntis_int.GetValue()
    def data_order(self): return self.preview.order_preview.order, self.preview.order_preview.tagfirst
    def tc_pairs(self): return self.tc_ch.checkbox.IsChecked()
    def labelling(self): return self.labelling_ch.GetSelection()
    def bolus_dur_type(self): return self.bolus_dur_ch.GetSelection()
    def repeats_type(self): return self.repeats_choice.GetSelection()
    def bolus_dur(self): 
        if self.bolus_dur_type() == 0: return [self.bolus_dur_num.GetValue(), ]
        else: return self.bolus_dur_list.GetValues()
    def repeats(self): return self.nrepeats_list.GetValues()
    def tis(self): 
        tis = self.ti_list.GetValues()
        if self.labelling() == 1:
            # For pASL TI = bolus_dur + PLD
            bolus_durs = self.bolus_dur()
            if len(bolus_durs) == 1: bolus_durs *= self.ntis()
            tis = [pld+bd for pld,bd in zip(tis, bolus_durs)]
        return tis
    def readout(self): return self.readout_ch.GetSelection()
    def time_per_slice(self): return self.time_per_slice_num.GetValue()
    def multiband(self): return self.multiband_cb.IsChecked()
    def slices_per_band(self): return self.slices_per_band_spin.GetValue()
    
    def set_default_dir(self, evt):
        """ 
        Bit of a hack - set the default dir for other file pickers to the same dir
        as the main data
        """
        f = self.data_picker.GetPath()
        try:
            self.nvols = nib.load(f).get_data().shape[3]
        except:
            # Failed to open file. Various reasons for this
            # but error will be displayed by the run handler
            self.nvols = -1

        d = os.path.dirname(f)
        for w in [self.structure.fsl_anat_picker,
                  self.structure.struc_image_picker,
                  self.structure.brain_image_picker,
                  self.calibration.calib_image_picker,
                  self.calibration.ref_tissue_mask_picker,
                  self.calibration.coil_image_picker,
                  self.analysis.mask_picker,
                  self.distcorr.calib_picker,
                  self.distcorr.fmap_picker,
                  self.distcorr.fmap_mag_picker,
                  self.distcorr.fmap_be_picker]:
            try:
                # WX version dependent...
                w.SetInitialDirectory(d)
            except:
                w.SetPath(d)
        self.update()

    def update(self, event=None):
        self.ti_list.set_size(self.ntis())
        self.bolus_dur_list.set_size(self.ntis())

        self.time_per_slice_num.Enable(self.readout() != 0)
        self.multiband_cb.Enable(self.readout() != 0)
        self.slices_per_band_spin.Enable(self.multiband() and self.readout() != 0)
        self.slices_per_band_label.Enable(self.multiband() and self.readout() != 0)

        self.bolus_dur_num.Enable(self.bolus_dur_type() == 0)
        self.bolus_dur_list.Enable(self.bolus_dur_type() == 1)

        if self.nvols == -1:
            # Couldn't find number of volumes in file
            self.nrepeats_list.set_size(self.ntis(), default=1)
        else:
            N = self.ntis() * (2 if self.tc_pairs() else 1)
            if self.repeats_type() == 0:
                # Clear out previous list of repeats to make sure number is fixed
                # Note run handler will detect error if not consistent with number of volumes
                self.nrepeats_list.set_size(0)
            self.nrepeats_list.set_size(self.ntis(), default=int(self.nvols / N))
            self.nrepeats_list.Enable(self.repeats_type() == 1)

        self.tc_ch.Enable(self.tc_pairs())
        self.update_groups()

        TabPage.update(self)

    def labelling_changed(self, event):
        if event.GetInt() == 0:
            self.bolus_dur_num.SetValue(0.7)
            self.ntis_int.label.SetLabel("Number of TIs")
            self.ti_list.label.SetLabel("TIs")
            self.preview.order_preview.tis_name="TIs"
            self.groups[0] = "TIs"
        else:
            self.bolus_dur_num.SetValue(1.8)
            self.ntis_int.label.SetLabel("Number of PLDs")
            self.ti_list.label.SetLabel("PLDs")
            self.preview.order_preview.tis_name="PLDs"
            self.groups[0] = "PLDs"
        self.analysis.labelling_changed(event.GetInt() == 0)
        self.update()

    def bolus_dur_changed(self, event):
        """ If constant bolus duration is changed, update the disabled list of
            bolus durations to match, to avoid any confusion """
        if self.bolus_dur_type() == 0:
            for c in range(self.bolus_dur_list.n): 
                self.bolus_dur_list.SetCellValue(0, c, str(self.bolus_dur()[0]))
        event.Skip()

    def update_groups(self, group1=True, group2=True):
        """
        This hideous code is to update the choices available in the ordering menus,
        hide the second if it is not relevant (1 PLD) and derive the data ordering
        string to pass to the data order preview
        """
        g1 = self.choice1.GetString(self.choice1.GetSelection())
        if self.choice2.IsShown():
            g2 = self.choice2.GetString(self.choice2.GetSelection())
        else:
            g2 = self.groups[0]
        if g1 == g2:
            g2 = self.groups[(self.groups.index(g1) + 1) % 3]

        self.choice2.Show()
        choices_st, choices_end = 0, 3

        # If data is not pairs, don't offer TC pairs options
        # Change selections to something else, but make sure they aren't the same
        if not self.tc_pairs():
            self.choice2.Hide() 
            if g1 == self.groups[2]: g1 = self.groups[0]
            if g2 in (g1, self.groups[2]): g2 = self.groups[1-self.groups.index(g1)]
            choices_end = 2
        
        # If only one TI/PLD, don't offer to group by TIs/PLDs
        # Change selections to something else. Second menu is set to group by PLDs
        # but is invisible now 
        if self.ntis() == 1:
            self.choice2.Hide()   
            if g1 == self.groups[0]: g1 = self.groups[1]
            g2 = self.groups[0]
            choices_st = 1
            
        group1_items = []
        group2_items = []
        for item in self.groups[choices_st:choices_end]:
            group1_items.append(item)
            if item != g1: group2_items.append(item)

        self.update_group_choice(self.choice1, group1_items, g1)
        self.update_group_choice(self.choice2, group2_items, g2)
        
        idx1 = self.groups.index(g1)
        idx2 = self.groups.index(g2)
        
        order = self.abbrevs[idx1]
        order += self.abbrevs[idx2]
        order += self.abbrevs[3-idx1-idx2]
        self.preview.order_preview.order = order

        # Need to do this as we may have unhidden the second menu
        self.GetSizer().Layout() 

    def update_group_choice(self, w, items, sel):
        w.Enable(False)
        w.Clear()
        w.AppendItems(items)
        w.SetSelection(w.FindString(sel))
        w.Enable(True)
