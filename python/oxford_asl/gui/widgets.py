"""
BASIL GUI for Oxford ASL - Base classes for pages in the tab notebook

Copyright (C) 2020 University of Oxford
"""
import os

import wx
import wx.grid

from . import OptionComponent

class TabPage(wx.Panel, OptionComponent):
    """
    Shared methods used by the various tab pages in the GUI
    """
    def __init__(self, app, notebook, title, tab_idx, num_tabs, name=None):
        OptionComponent.__init__(self, app)
        wx.Panel.__init__(self, parent=notebook, id=wx.ID_ANY)
        self.notebook = notebook
        self.tab_idx = tab_idx
        self.num_tabs = num_tabs
        self.sizer = wx.GridBagSizer(vgap=5, hgap=5)
        self.row = 0
        self.title = title
        if name is None:
            self.name = title.lower()
        else:
            self.name = name

    def add_next_prev_btn(self):
        """
        Add next/previous buttons
        """
        if self.tab_idx < self.num_tabs-1:
            next_btn = wx.Button(self, label="Next", id=wx.ID_FORWARD)
            next_btn.Bind(wx.EVT_BUTTON, self._next)
        else:
            next_btn = wx.StaticText(self, label="")

        if self.tab_idx > 0:
            prev_btn = wx.Button(self, label="Previous", id=wx.ID_BACKWARD)
            prev_btn.Bind(wx.EVT_BUTTON, self._prev)
        else:
            prev_btn = wx.StaticText(self, label="")

        self.pack(" ")
        self.sizer.AddGrowableRow(self.row-1, 1)
        self.sizer.Add(prev_btn, pos=(self.row, 0), border=10, flag=wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_LEFT)
        self.sizer.Add(wx.StaticText(self, label=""), pos=(self.row, 1), border=10, flag=wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_LEFT)
        self.sizer.Add(wx.StaticText(self, label=""), pos=(self.row, 2), border=10, flag=wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_LEFT)
        self.sizer.Add(next_btn, pos=(self.row, 3), border=10, flag=wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_RIGHT)

    def _next(self, _evt):
        self.notebook.SetSelection(self.tab_idx+1)

    def _prev(self, _evt):
        self.notebook.SetSelection(self.tab_idx-1)

    def pack(self, label, *widgets, **kwargs):
        """
        Add a horizontal line to the tab with a label and series of widgets

        If label is empty, first widget is used instead (usually to provide a checkbox)
        """
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

        for widget in widgets:
            span = (1, 1)
            widget.label = text
            if hasattr(widget, "span"):
                span = (1, widget.span)
            widget.SetFont(font)
            widget.Enable(col == 0 or kwargs.get("enable", True))
            self.sizer.Add(widget, pos=(self.row, col), border=border, flag=wx.ALIGN_CENTRE_VERTICAL | wx.EXPAND | wx.LEFT, span=span)
            col += span[1]
        self.row += 1

    def file_picker(self, label, pick_dir=False, handler=None, optional=False, initial_on=False, pack=True, **kwargs):
        """
        Add a file picker to the tab
        """
        if not handler:
            handler = self.state_changed
        if pick_dir:
            picker = wx.DirPickerCtrl(self, style=wx.DIRP_USE_TEXTCTRL)
            picker.Bind(wx.EVT_DIRPICKER_CHANGED, handler)
        else:
            picker = wx.FilePickerCtrl(self)
            picker.Bind(wx.EVT_FILEPICKER_CHANGED, handler)
        picker.span = 2
        if optional:
            checkbox = wx.CheckBox(self, label=label)
            checkbox.SetValue(initial_on)
            checkbox.Bind(wx.EVT_CHECKBOX, self._checkbox_toggle_cb(checkbox, picker))
            picker.checkbox = checkbox
            if pack:
                self.pack("", checkbox, picker, enable=initial_on, **kwargs)
        elif pack:
            self.pack(label, picker, **kwargs)

        return picker

    def set_picker_dir(self, filepicker, fpath):
        """
        Set the initial directory for a file or dir picker

        :param fpath: Path to file or directory
        """
        dirname = os.path.dirname(fpath)
        try:
            filepicker.SetInitialDirectory(dirname)
        except AttributeError:
            # WX version dependent - so try alternate name
            filepicker.SetPath(dirname)

    def _checkbox_toggle_cb(self, checkbox, widget):
        def _toggled(_event):
            widget.Enable(checkbox.IsChecked())
            self.state_changed(_event)
        return _toggled

    def choice(self, label, choices, initial=0, optional=False, initial_on=False, handler=None, pack=True, **kwargs):
        """
        Add a widget to choose from a fixed set of options
        """
        if not handler:
            handler = self.state_changed
        choice = wx.Choice(self, choices=choices)
        choice.SetSelection(initial)
        choice.Bind(wx.EVT_CHOICE, handler)
        if optional:
            checkbox = wx.CheckBox(self, label=label)
            checkbox.SetValue(initial_on)
            checkbox.Bind(wx.EVT_CHECKBOX, self.state_changed)
            choice.checkbox = checkbox
            if pack:
                self.pack("", checkbox, choice, enable=initial_on, **kwargs)
        elif pack:
            self.pack(label, choice, **kwargs)
        return choice

    def number(self, label, handler=None, **kwargs):
        """
        Add a widget to choose a floating point number
        """
        if not handler:
            handler = self.state_changed
        num = NumberChooser(self, changed_handler=handler, **kwargs)
        num.span = 2
        self.pack(label, num, **kwargs)
        return num

    def integer(self, label, handler=None, pack=True, **kwargs):
        """
        Add a widget to choose an integer
        """
        if not handler:
            handler = self.state_changed
        spin = wx.SpinCtrl(self, **kwargs)
        spin.SetValue(kwargs.get("initial", 0))
        spin.Bind(wx.EVT_SPINCTRL, handler)
        if pack:
            self.pack(label, spin)
        return spin

    def checkbox(self, label, initial=False, handler=None, **kwargs):
        """
        Add a simple on/off option
        """
        checkbox = wx.CheckBox(self, label=label)
        checkbox.span = 2
        checkbox.SetValue(initial)
        if handler:
            checkbox.Bind(wx.EVT_CHECKBOX, handler)
        else:
            checkbox.Bind(wx.EVT_CHECKBOX, self.state_changed)
        self.pack("", checkbox, **kwargs)
        return checkbox

    def section(self, label):
        """
        Add a section heading
        """
        self.pack(label, bold=True)

class WhitePaperCompatibility(OptionComponent, wx.Panel):
    """
    Widget for displaying White Paper compatibility and resetting defaults

    White Paper compatibility is defined by:
      - Voxelwise calibration
      - Fixed T1 of 1.65
      - Arterial transit time of 0

    In addition, white paper mode was defined for single-PLD data and therefore no inference of ATT
    """

    def __init__(self, app, parent):
        OptionComponent.__init__(self, app)
        wx.Panel.__init__(self, parent)
        self._wp = False
        self._issues = []
        self._img_compatible = os.path.join(os.path.abspath(os.path.dirname(__file__)), "wp_tick.png")
        self._img_incompatible = os.path.join(os.path.abspath(os.path.dirname(__file__)), "wp_cross.png")
        
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self._img = wx.Panel(self, size=(48, 48))
        self.hbox.Add(self._img)

        status_panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)
        status_panel.SetSizer(vbox)

        self._label = wx.StaticText(status_panel, label="")
        vbox.Add(self._label)

        btn_panel = wx.Panel(status_panel)
        btn_hbox = wx.BoxSizer(wx.HORIZONTAL)
        self._reset_btn = wx.Button(btn_panel, label="Make compatible")
        self._reset_btn.Bind(wx.EVT_BUTTON, self._make_compatible)
        btn_hbox.Add(self._reset_btn)
        self._view_btn = wx.Button(btn_panel, label="View issues")
        self._view_btn.Bind(wx.EVT_BUTTON, self._view_issues)
        btn_hbox.Add(self._view_btn)
        btn_panel.SetSizer(btn_hbox)
        vbox.Add(btn_panel)
        
        self._warn = wx.StaticText(status_panel, label="")
        self._warn.SetForegroundColour(wx.Colour(255, 0, 0))
        vbox.Add(self._warn)

        self.hbox.Add(status_panel, wx.EXPAND)
        self.SetSizer(self.hbox)

    def options(self):
        return {"wp" : self._wp}

    def option_changed(self, options, key, value):
        self._issues = self._get_issues(options)

        if len(self._issues) == 0:
            wx.StaticBitmap(self._img, -1, wx.Bitmap(self._img_compatible, wx.BITMAP_TYPE_ANY))
            self._label.SetLabel("Analysis is compatible with white paper (Alsop et al 2004)")
            self._wp = True
        else:
            wx.StaticBitmap(self._img, -1, wx.Bitmap(self._img_incompatible, wx.BITMAP_TYPE_ANY))
            self._label.SetLabel("Analysis is not compatible with white paper (Alsop et al 2004)")
            self._wp = False

        if not self._issues and options["_ntis"] > 1:
            self._warn.SetLabel("Warning: White paper was intended single TI/PLD data")
        else:
            self._warn.SetLabel("")

        self._reset_btn.Enable(not self._wp)
        self._view_btn.Enable(not self._wp)

    def _get_issues(self, options):
        issues = []
        if "c" not in options:
            issues.append(("Calibration image is not provided - results will not be quantitative", "Calibration will be turned on"))
        elif options["cmethod"] != "voxel":
            issues.append(("White paper specifies voxelwise calibration - you are using reference region", "Calibration method will be switched to voxelwise"))

        return issues
            
    def _make_compatible(self, _event=None):
        fix_list = "\n".join(["%i. %s" % (idx+1, issue[1]) for idx, issue in enumerate(self._issues)])
        message = "The following changes will be made:\n\n" + fix_list
        ret = wx.MessageBox(message, "Make analysis WP compatible", wx.OK | wx.CANCEL | wx.ICON_INFORMATION)  
        if ret == wx.OK:
            self._wp = True
            self.state_changed()

    def _view_issues(self, _event=None):
        issue_list = "\n".join(["%i. %s" % (idx+1, issue[0]) for idx, issue in enumerate(self._issues)])
        message = "The current analysis is not WP compatible for the following reasons:\n\n" + issue_list
        wx.MessageBox(message, "Reasons for WP incompatibility", wx.OK | wx.ICON_INFORMATION)  

class NumberChooser(wx.Panel):
    """
    Widget for choosing a floating point number
    """

    def __init__(self, parent, label=None, minval=0, maxval=1, initial=0.5, step=0.1, digits=2, changed_handler=None):
        super(NumberChooser, self).__init__(parent)
        self.minval, self.orig_minval, self.maxval, self.orig_maxval = minval, minval, maxval, maxval
        self.handler = changed_handler
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        if label is not None:
            self.label = wx.StaticText(self, label=label)
            self.hbox.Add(self.label, proportion=0, flag=wx.ALIGN_CENTRE_VERTICAL)
        # Set a very large maximum as we want to let the user override the default range
        #self.spin = wx.SpinCtrl(self, minval=0, maxval=100000, initial=initial)
        #self.spin.Bind(wx.EVT_SPINCTRL, self._spin_changed)
        self.spin = wx.SpinCtrlDouble(self, min=0, max=100000, inc=step, initial=initial)
        self.spin.SetDigits(digits)
        self.spin.Bind(wx.EVT_SPINCTRLDOUBLE, self._spin_changed)
        self.slider = wx.Slider(self, value=initial, minValue=0, maxValue=100)
        self.slider.SetValue(100*(initial-self.minval)/(self.maxval-self.minval))
        self.slider.Bind(wx.EVT_SLIDER, self._slider_changed)
        self.hbox.Add(self.slider, proportion=1, flag=wx.EXPAND | wx.ALIGN_CENTRE_VERTICAL)
        self.hbox.Add(self.spin, proportion=0, flag=wx.EXPAND | wx.ALIGN_CENTRE_VERTICAL)
        self.SetSizer(self.hbox)

    def GetValue(self):
        """
        Get the currently selected number
        """
        return self.spin.GetValue()

    def SetValue(self, val):
        """
        Set the selected number
        """
        self.spin.SetValue(val)
        self.slider.SetValue(100*(val-self.minval)/(self.maxval-self.minval))

    def _slider_changed(self, event):
        slider_pos = event.GetInt()
        val = self.minval + (self.maxval-self.minval)*float(slider_pos)/100
        self.spin.SetValue(val)
        if self.handler:
            self.handler(event)
        event.Skip()

    def _spin_changed(self, event):
        """ If user sets the spin outside the current range, update the slider range
        to match. However if they go back inside the current range, revert to this for
        the slider"""
        val = event.GetValue()
        if val < self.minval:
            self.minval = val
        elif val > self.orig_minval:
            self.minval = self.orig_minval
        if val > self.maxval:
            self.maxval = val
        elif val < self.orig_maxval:
            self.maxval = self.maxval
        self.slider.SetValue(100*(val-self.minval)/(self.maxval-self.minval))
        if self.handler:
            self.handler()
        event.Skip()

class NumberList(wx.grid.Grid):
    """
    Widget for specifying a list of numbers
    """

    def __init__(self, parent, size, default=1.8):
        super(NumberList, self).__init__(parent, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, 0)
        self.size = 0
        self.default = default
        self.CreateGrid(1, 0)
        self.SetRowLabelSize(0)
        self.SetColLabelSize(0)
        self.SetSize(size)
        self.Bind(wx.EVT_SIZE, self._on_size)

    def GetValues(self):
        """
        :return: Sequence of values in the list
        """
        try:
            return [float(self.GetCellValue(0, c)) for c in range(self.size)]
        except ValueError:
            raise RuntimeError("Non-numeric values in number list")

    def SetSize(self, size, default=None):
        """
        Set the size of the number list

        :param size: Number of items in list
        :param default: Default value to use for newly created columns. If not specified
                        uses default defined in constructor
        """
        if default is None:
            if self.size == 0:
                default = self.default
            else:
                default = self.GetCellValue(0, self.size-1)
        if size > self.size:
            self.AppendCols(size - self.size)
            for col in range(self.size, size):
                self.SetCellValue(0, col, str(default))
        elif size < self.size:
            self.DeleteCols(size, self.size-size)
        self.size = size
        self._resize_cols()

    def _resize_cols(self):
        if self.size == 0:
            return

        width, _height = self.GetClientSize()
        col_width = width / self.size
        for i in range(self.size):
            self.SetColSize(i, col_width)

    def _on_size(self, event):
        event.Skip()
        self._resize_cols()
