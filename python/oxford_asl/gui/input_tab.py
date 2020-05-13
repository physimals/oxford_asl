"""
BASIL GUI for Oxford ASL - Input data/acquisition options

Copyright (C) 2020 University of Oxford
"""
import wx
import wx.grid

from . import get_nvols, OptionError
from .widgets import TabPage, NumberChooser, NumberList

class AslInputOptions(TabPage):
    """
    Tab page containing input data options
    """

    IAF_CHOICES = ["tc", "ct", "diff"]

    def __init__(self, app, parent, idx, n):
        TabPage.__init__(self, app, parent, "Input Data", idx, n, name="input")
        self._nvols = -1 # Cached to improve responsiveness

        self.section("Data contents")
        self.data_picker = self.file_picker("Input Image")
        self.ntis_int = self.integer("Number of PLDs", min=1, max=100, initial=1)
        self.repeats_choice = self.choice("Repeats", choices=["Fixed", "Variable"])

        self.section("Data order")
        self.ibf_choice = self.choice("Volumes grouped by", choices=["Repeats", "TIs/PLDs"])
        self.iaf_choice = self.choice("Label/Control pairing", choices=["Label then control", "Control then label", "Pre-subtracted"])

        self.section("Acquisition parameters")
        self.labelling_ch = self.choice("Labelling", choices=["pASL", "cASL/pcASL"], initial=1)
        self.bolus_dur_ch = wx.Choice(self, choices=["Constant", "Variable"])
        self.bolus_dur_ch.SetSelection(0)
        self.bolus_dur_ch.Bind(wx.EVT_CHOICE, self.state_changed)
        self.bolus_dur_num = NumberChooser(self, minval=0, maxval=2.5, step=0.1, initial=1.8, changed_handler=self.state_changed)
        self.bolus_dur_num.span = 2
        self.pack("Bolus duration (s)", self.bolus_dur_ch, self.bolus_dur_num)

        self.bolus_dur_list = NumberList(self, 1)
        self.bolus_dur_list.span = 3
        self.bolus_dur_list.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self.state_changed)
        self.pack("Bolus durations (s)", self.bolus_dur_list, enable=False)

        self.ti_list = NumberList(self, 1)
        self.ti_list.span = 3
        self.ti_list.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self.state_changed)
        self.pack("PLDs (s)", self.ti_list)

        self.nrepeats_list = NumberList(self, 1, default=1)
        self.nrepeats_list.span = 3
        self.nrepeats_list.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self.state_changed)
        self.pack("Repeats", self.nrepeats_list, enable=False)

        self.readout_ch = wx.Choice(self, choices=["3D (eg GRASE)", "2D multi-slice (eg EPI)"])
        self.readout_ch.SetSelection(0)
        self.readout_ch.Bind(wx.EVT_CHOICE, self.state_changed)
        self.time_per_slice_num = NumberChooser(self, label="Time per slice (ms)", minval=0, maxval=50, step=1, initial=10.0, changed_handler=self.state_changed)
        self.time_per_slice_num.span = 2
        self.pack("Readout", self.readout_ch, self.time_per_slice_num)
        self.time_per_slice_num.Enable(False)

        self.multiband_cb = wx.CheckBox(self, label="Multi-band")
        self.multiband_cb.Bind(wx.EVT_CHECKBOX, self.state_changed)
        self.slices_per_band_spin = wx.SpinCtrl(self, min=1, max=100, initial=5)
        self.slices_per_band_spin.Bind(wx.EVT_SPINCTRL, self.state_changed)
        self.slices_per_band_label = wx.StaticText(self, label="slices per band")
        self.pack("", self.multiband_cb, self.slices_per_band_spin, self.slices_per_band_label, enable=False)
        self.multiband_cb.Enable(False)

        self.sizer.AddGrowableCol(2, 1)
        self.SetSizer(self.sizer)
        self.add_next_prev_btn()
        self.Layout()

    def options(self):
        options = {
            "i"          : self.data_picker.GetPath(),
            "iaf"        : self.IAF_CHOICES[self.iaf_choice.GetSelection()],
            "ibf"        : "rpt" if self.ibf_choice.GetSelection() == 0 else "tis",
            "casl"       : bool(self.labelling_ch.GetSelection()),
            "bolus"      : self.bolus_dur_list.GetValues() if bool(self.bolus_dur_ch.GetSelection()) else self.bolus_dur_num.GetValue(),
            "rpts"       : [int(v) for v in self.nrepeats_list.GetValues()],
            "slicedt"    : self.time_per_slice_num.GetValue() / 1000 if bool(self.readout_ch.GetSelection()) else None,
            "sliceband"  : self.slices_per_band_spin.GetValue() if self.multiband_cb.IsChecked() else None,
            "_ntis"      : self.ntis_int.GetValue(),
            "_nvols"     : self._nvols,
            "_ntc"       : 1 if self.iaf_choice.GetSelection() == 2 else 2,
            "_fixbolus"  : not bool(self.bolus_dur_ch.GetSelection()),
            "_fixrpts"   : not bool(self.repeats_choice.GetSelection())
        }
        options["tis"] = self._calculate_tis(options)
        return options

    def check_options(self, options):
        self._check_image("Input data", options["i"], can_be_none=False)
            
        if min(options["rpts"]) < 0:
            raise OptionError("Invalid value in repeats list")
        if not isinstance(options["bolus"], float) and min(options["bolus"]) < 0:
            raise OptionError("Invalid bolus duration")
        if min(options["tis"]) < 0:
            raise OptionError("Invalid value in TIs/PLDs list")
        
        expected_vols = options["_ntc"] * sum(options["rpts"])
        if expected_vols != options["_nvols"] and options["_nvols"] > 0:
            if options["_fixrpts"]:
                raise OptionError("%i volumes in data - inconsistent with number of TIs/PLDs" % options["_nvols"])
            else:
                raise OptionError("%i volumes in data - inconsistent with number of TIs/PLDs and total number of repeats" % options["_nvols"])

    def option_changed(self, options, key, value):
        if key == "i":
            self._nvols = get_nvols(self.data_picker.GetPath())
        if key == "_ntis":
            self.ti_list.SetNumValues(value)
            self.bolus_dur_list.SetNumValues(value)
            self.nrepeats_list.SetNumValues(value)
            self._update_fixed_repeats(options)

        if key in ("_nvols", "_ntc"):
            self._update_fixed_repeats(options)

        if key == "_fixbolus":
            self.bolus_dur_num.Enable(value)
            self.bolus_dur_list.Enable(not value)

        if key == "_fixrpts":
            self.nrepeats_list.Enable(not value)

        if key == "casl":
            self._labelling_changed(options)

        if key == "bolus":
            self._bolus_changed(options)

        if key in ("slicedt", "sliceband"):
            self.time_per_slice_num.Enable(options["slicedt"] is not None)
            self.multiband_cb.Enable(options["slicedt"] is not None)
            self.slices_per_band_spin.Enable(options["sliceband"] is not None and options["slicedt"] is not None)
            self.slices_per_band_label.Enable(options["sliceband"] is not None and options["slicedt"] is not None)

    def _calculate_tis(self, options):
        tis = self.ti_list.GetValues()
        if options["casl"]:
            # For pASL TI = bolus_dur + PLD
            if options["_fixbolus"]:
                bolus_durs = [options["bolus"]] * len(tis)
            else:
                bolus_durs = options["bolus"]

            tis = [pld+bolus_dur for pld, bolus_dur in zip(tis, bolus_durs)]
        return tis

    def _update_fixed_repeats(self, options):
        if options["_fixrpts"] and options["_nvols"] > 0:
            # Clear out previous list of repeats and reset to make sure number is fixed
            # Note check_options will detect error if not consistent with number of volumes
            vols_per_repeat = options["_ntis"] * options["_ntc"]
            self.nrepeats_list.SetNumValues(0)
            self.nrepeats_list.SetNumValues(options["_ntis"], default=int(options["_nvols"] / vols_per_repeat))

    def _labelling_changed(self, options):
        if options["casl"]:
            self.bolus_dur_num.SetValue(1.8)
            times_name = "PLDs"
        else:
            self.bolus_dur_num.SetValue(0.7)
            times_name = "TIs"

        self.ntis_int.label.SetLabel("Number of %s" % times_name)
        self.ti_list.label.SetLabel(times_name)

    def _bolus_changed(self, options):
        # If constant bolus duration is changed, update the disabled list of
        # bolus durations to match, to avoid any confusion
        if options["_fixbolus"]:
            for col in range(self.bolus_dur_list.size):
                self.bolus_dur_list.SetCellValue(0, col, str(options["bolus"]))
