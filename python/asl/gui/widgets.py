import sys
import os
import colorsys

import wx
import wx.grid

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure

import numpy as np

class TabPage(wx.Panel):
    """
    Shared methods used by the various tab pages in the GUI
    """
    def __init__(self, parent, title, idx, n, name=None):
        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)
        self.notebook = parent
        self.idx = idx
        self.n = n
        self.sizer = wx.GridBagSizer(vgap=5, hgap=5)
        self.row = 0
        self.title = title
        if name is None:
            self.name = title.lower()
        else:
            self.name = name
            
    def next_prev(self):
        """
        Add next/previous buttons
        """
        if self.idx < self.n-1: 
            self.next_btn = wx.Button(self, label="Next", id=wx.ID_FORWARD)
            self.next_btn.Bind(wx.EVT_BUTTON, self.next)
        else:
            self.next_btn = wx.StaticText(self, label="")

        if self.idx > 0:
            self.prev_btn = wx.Button(self, label="Previous", id=wx.ID_BACKWARD)
            self.prev_btn.Bind(wx.EVT_BUTTON, self.prev)
        else:
            self.prev_btn = wx.StaticText(self, label="")

        self.pack(" ")
        self.sizer.AddGrowableRow(self.row-1, 1)
        self.sizer.Add(self.prev_btn, pos=(self.row, 0), border=10, flag=wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_LEFT)
        self.sizer.Add(wx.StaticText(self, label=""), pos=(self.row, 1), border=10, flag=wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_LEFT)
        self.sizer.Add(wx.StaticText(self, label=""), pos=(self.row, 2), border=10, flag=wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_LEFT)
        self.sizer.Add(self.next_btn, pos=(self.row, 3), border=10, flag=wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_RIGHT)

    def next(self, evt):
        self.notebook.SetSelection(self.idx+1)

    def prev(self, evt):
        self.notebook.SetSelection(self.idx-1)

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
        """
        Add a file picker to the tab
        """
        if not handler: handler = self.update
        if dir: 
            picker = wx.DirPickerCtrl(self, style=wx.DIRP_USE_TEXTCTRL)
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
        """
        Add a widget to choose from a fixed set of options
        """
        if not handler: handler = self.update
        ch = wx.Choice(self, choices=choices)
        ch.SetSelection(initial)
        ch.Bind(wx.EVT_CHOICE, handler)
        if optional:
            cb = wx.CheckBox(self, label=label)
            cb.SetValue(initial_on)
            cb.Bind(wx.EVT_CHECKBOX, self.update)
            ch.checkbox = cb
            if pack: self.pack("", cb, ch, enable=initial_on, **kwargs)
        elif pack:
            self.pack(label, ch, **kwargs)
        return ch

    def number(self, label, handler=None, **kwargs):
        """
        Add a widget to choose a floating point number
        """
        if not handler: handler = self.update
        num = NumberChooser(self, changed_handler=handler, **kwargs)
        num.span = 2
        self.pack(label, num, **kwargs)
        return num

    def integer(self, label, handler=None, pack=True, **kwargs):
        """
        Add a widget to choose an integer
        """
        if not handler: handler = self.update
        spin = wx.SpinCtrl(self, **kwargs)
        spin.SetValue(kwargs.get("initial", 0))
        spin.Bind(wx.EVT_SPINCTRL, handler)
        if pack: self.pack(label, spin)
        return spin

    def checkbox(self, label, initial=False, handler=None, **kwargs):
        """
        Add a simple on/off option
        """
        cb = wx.CheckBox(self, label=label)
        cb.span=2
        cb.SetValue(initial)
        if handler: cb.Bind(wx.EVT_CHECKBOX, handler)
        else: cb.Bind(wx.EVT_CHECKBOX, self.update)
        self.pack("", cb, **kwargs)
        return cb

    def section(self, label):
        """
        Add a section heading
        """
        self.pack(label, bold=True)

    def update(self, evt=None):
        """
        Update the run module, i.e. when options have changed
        """
        if hasattr(self, "run"): 
            self.run.update()
            if hasattr(self, "preview"): self.preview.run = self.run

class NumberChooser(wx.Panel):
    """
    Widget for choosing a floating point number
    """

    def __init__(self, parent, label=None, min=0, max=1, initial=0.5, step=0.1, digits=2, changed_handler=None):
        super(NumberChooser, self).__init__(parent)
        self.min, self.orig_min, self.max, self.orig_max = min, min, max, max
        self.handler = changed_handler
        self.hbox=wx.BoxSizer(wx.HORIZONTAL)
        if label is not None:
            self.label = wx.StaticText(self, label=label)
            self.hbox.Add(self.label, proportion=0, flag=wx.ALIGN_CENTRE_VERTICAL)
        # Set a very large maximum as we want to let the user override the default range
        #self.spin = wx.SpinCtrl(self, min=0, max=100000, initial=initial)
        #self.spin.Bind(wx.EVT_SPINCTRL, self.spin_changed)
        self.spin = wx.SpinCtrlDouble(self, min=0, max=100000, inc=step, initial=initial)
        self.spin.SetDigits(digits)
        self.spin.Bind(wx.EVT_SPINCTRLDOUBLE, self.spin_changed)
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
        event.Skip()

    def spin_changed(self, event):
        """ If user sets the spin outside the current range, update the slider range
        to match. However if they go back inside the current range, revert to this for
        the slider"""
        val = event.GetValue()
        if val < self.min: 
            self.min = val
        elif val > self.orig_min:
            self.min = self.orig_min
        if val > self.max: 
            self.max = val
        elif val < self.orig_max:
            self.max = self.orig_max
        self.slider.SetValue(100*(val-self.min)/(self.max-self.min))
        if self.handler: self.handler(event)
        event.Skip()
        
class NumberList(wx.grid.Grid):
    """
    Widget for specifying a list of numbers
    """

    def __init__(self, parent, n, default=1.8):
        super(NumberList, self).__init__(parent, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, 0 )
        self.n=0
        self.default = default
        self.CreateGrid(1, 0)
        self.SetRowLabelSize(0)
        self.SetColLabelSize(0)
        self.set_size(n)
        self.Bind(wx.EVT_SIZE, self.on_size)

    def GetValues(self):
        try:
            return [float(self.GetCellValue(0, c)) for c in range(self.n)]
        except ValueError as e:
            raise RuntimeError("Non-numeric values in number list")
            
    def set_size(self, n):
        if self.n == 0: default = self.default
        else: default = self.GetCellValue(0, self.n-1)
        if n > self.n:
            self.AppendCols(n - self.n)
            for c in range(self.n, n): self.SetCellValue(0, c, str(default))
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

class PreviewPanel(wx.Panel):
    """
    Panel providing a simple image preview for the output of ASL_FILE.

    Used so user can check their choice of data grouping/ordering looks right
    """
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, size=wx.Size(300, 600))
        self.data = None
        self.run = None
        self.slice = -1
        self.nslices = 1
        self.view = 0
        self.figure = Figure(figsize=(3.5, 3.5), dpi=100, facecolor='black')
        self.axes = self.figure.add_subplot(111, facecolor='black')
        self.axes.get_xaxis().set_ticklabels([])
        self.axes.get_yaxis().set_ticklabels([])          
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.canvas.mpl_connect('scroll_event', self.scroll)
        self.canvas.mpl_connect('button_press_event', self.view_change)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        font = self.GetFont()
        font.SetWeight(wx.BOLD)
        text = wx.StaticText(self, label="Data preview - perfusion weighted image")
        text.SetFont(font)
        self.sizer.AddSpacer(10)
        self.sizer.Add(text, 0)   
        self.sizer.Add(self.canvas, 2, border=5, flag = wx.EXPAND | wx.ALL)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(wx.StaticText(self, label="Use scroll wheel to change slice, double click to change view"), 0, flag=wx.ALIGN_CENTRE_VERTICAL)      
        self.update_btn = wx.Button(self, label="Update")
        self.update_btn.Bind(wx.EVT_BUTTON, self.update)
        hbox.Add(self.update_btn)
        self.sizer.Add(hbox)

        self.sizer.AddSpacer(10)
        text = wx.StaticText(self, label="Data order preview")
        text.SetFont(font)
        self.sizer.Add(text, 0)
        self.order_preview = AslDataPreview(self, 1, 1, True, "trp", True)
        self.sizer.Add(self.order_preview, 2, wx.EXPAND)
        self.SetSizer(self.sizer)
        self.Layout()

    def update(self, evt):
        """
        Update the preview. This is called explicitly when the user clicks the update
        button as it involves calling ASL_FILE and may be slow
        """
        self.data = None
        if self.run is not None:
            self.data = self.run.get_preview_data()
            # If multi-TI data, take mean over volumes
            if self.data is not None and len(self.data.shape) == 4:
                self.data = np.mean(self.data, axis=3)
                
        if self.data is not None:
            self.view = 0
            self.init_view()
        self.redraw()

    def init_view(self):
        self.nslices = self.data.shape[2-self.view]
        self.slice = int(self.nslices / 2)
        self.redraw()

    def redraw(self):
        """
        Redraw the preview image
        """
        self.axes.clear() 
        if self.data is None: return
        if self.view == 0:
            sl = self.data[:,:,self.slice]
        elif self.view == 1:
            sl = self.data[:,self.slice,:]
        else:
            sl = self.data[self.slice,:,:]

        i = self.axes.imshow(sl.T, interpolation="nearest", vmin=sl.min(), vmax=sl.max())
        self.axes.set_ylim(self.axes.get_ylim()[::-1])
        i.set_cmap("gray")
        self.Layout()

    def view_change(self, event):
        """
        Called on mouse click event. Double click changes the view direction and redraws
        """
        if self.data is None: return
        if event.dblclick:
            self.view = (self.view + 1) % 3
            self.init_view()
            self.redraw()

    def scroll(self, event):
        """
        Called on mouse scroll wheel to move through the slices in the current view
        """
        if event.button == "up":
            if self.slice != self.nslices-1: self.slice += 1
        else:
            if self.slice != 0: self.slice -= 1
        self.redraw()
            
class AslDataPreview(wx.Panel):
    """
    Widget to display a preview of the data ordering selected (i.e. how the volumes in a 4D
    dataset map to TIs, repeats and tag/control pairs)
    """

    def __init__(self, parent, n_tis, n_repeats, tc_pairs, order, tagfirst):
        wx.Panel.__init__(self, parent, size=wx.Size(300, 300))
        self.SetBackgroundStyle(wx.BG_STYLE_CUSTOM)
        self.Bind(wx.EVT_SIZE, self.on_size)
        self.Bind(wx.EVT_PAINT, self.on_paint)
        self.n_tis = n_tis
        self.n_repeats = int(n_repeats)
        self.tc_pairs = tc_pairs
        self.tagfirst = tagfirst
        self.order = order
        self.tis_name = "PLDs"

    def on_size(self, event):
        event.Skip()
        self.Refresh()

    def get_col(self, pos, ti):
        if ti: h = 170.0/255
        else: 
            h = 90.0/255
        s, v = 0.5, 0.95 - pos/2
        r,g,b = colorsys.hsv_to_rgb(h, s, v)
        return wx.Colour(int(r*255), int(g*255), int(b*255))

    def on_paint(self, event):
        w, h = self.GetClientSize()
        N = self.n_tis * int(self.n_repeats)
        if self.tc_pairs: N *= 2
        dc = wx.AutoBufferedPaintDC(self)
        dc.Clear()

        leg_width = (w-100)/4
        leg_start = 50

        dc.SetBrush(wx.TRANSPARENT_BRUSH)
        rect = wx.Rect(leg_start, 20, leg_width/4, 20)
        dc.GradientFillLinear(rect, self.get_col(0, True), self.get_col(1.0, True), wx.EAST)
        dc.DrawRectangle(*rect.Get())
        dc.DrawText(self.tis_name, leg_start+leg_width/3, 20)

        rect = wx.Rect(leg_start+leg_width, 20, leg_width/4, 20)
        dc.GradientFillLinear(rect, self.get_col(0, False), self.get_col(1.0, False), wx.EAST)
        dc.DrawRectangle(*rect.Get())
        dc.DrawText("Repeats", leg_start+4*leg_width/3, 20)

        if self.tc_pairs:
            dc.SetBrush(wx.TRANSPARENT_BRUSH)
            dc.DrawRectangle(leg_start+leg_width*2, 20, leg_width/4, 20)
            dc.DrawText("Label", leg_start+7*leg_width/3, 20)

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
                    seq += [i+2,] * (int(self.n_repeats) - 1)
        
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
