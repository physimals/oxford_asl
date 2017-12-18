import wx
import wx.grid

from .widgets import TabPage

class StructureTab(TabPage):
    EXISTING_FSLANAT = 0
    NEW_FSLANAT = 1
    INDEP_STRUC = 2
    NONE = 3

    TRANS_MATRIX = 0
    TRANS_IMAGE = 1
    TRANS_FSLANAT = 2

    """
    Tab page containing options for structural space transformation
    """
    def __init__(self, parent, idx, n):
        TabPage.__init__(self, parent, "Structure", idx, n, name="structure")

        self.section("Structure")

        self.struc_ch = wx.Choice(self, choices=["Existing FSL_ANAT output", "Run FSL_ANAT on structural image", "Independent structural data", "None"])
        self.struc_ch.SetSelection(self.NONE)
        self.struc_ch.Bind(wx.EVT_CHOICE, self.update)
        self.struc_ch.span = 2
        self.pack("Structural data from", self.struc_ch)
 
        self.fsl_anat_picker = self.file_picker("Existing FSL_ANAT directory", dir=True)
        self.struc_image_picker = self.file_picker("Structural Image")
        self.brain_image_picker = self.file_picker("Brain image", optional=True)
         
        self.section("Registration")

        self.transform_choices = ["Use matrix", "Use warp image", "Use FSL_ANAT"]

        self.transform_cb = wx.CheckBox(self, label="Transform to standard space")
        self.transform_cb.Bind(wx.EVT_CHECKBOX, self.update)
        self.transform_ch = wx.Choice(self, choices=self.transform_choices)
        self.transform_ch.SetSelection(self.TRANS_FSLANAT)
        self.transform_ch.Bind(wx.EVT_CHOICE, self.update)
        self.transform_picker = wx.FilePickerCtrl(self)
        self.transform_picker.Bind(wx.EVT_FILEPICKER_CHANGED, self.update)
        self.pack("", self.transform_cb, self.transform_ch, self.transform_picker, enable=False)

        self.sizer.AddGrowableCol(2, 1)
        self.SetSizer(self.sizer)
        self.next_prev()

    def existing_fsl_anat(self): 
        if self.struc_ch.GetSelection() == self.EXISTING_FSLANAT:
            return self.fsl_anat_picker.GetPath()
        else: 
            return None

    def run_fsl_anat(self):
        return self.struc_ch.GetSelection() == self.NEW_FSLANAT
        
    def struc_image(self): 
        if self.struc_ch.GetSelection() in (self.NEW_FSLANAT, self.INDEP_STRUC):
            return self.struc_image_picker.GetPath()
        else: 
            return None
            
    def struc_image_brain(self): 
        if self.struc_ch.GetSelection() == self.INDEP_STRUC and self.brain_image_picker.checkbox.IsChecked():
            return self.brain_image_picker.GetPath()
        else: return None

    def transform(self): 
        return self.transform_cb.IsChecked()

    def transform_type(self): 
        return self.transform_ch.GetSelection()

    def transform_file(self): 
        return self.transform_picker.GetPath()
    
    def update(self, event=None):
        mode = self.struc_ch.GetSelection()
        self.fsl_anat_picker.Enable(mode == self.EXISTING_FSLANAT)
        self.struc_image_picker.Enable(mode in (self.NEW_FSLANAT, self.INDEP_STRUC))

        self.brain_image_picker.checkbox.Enable(mode == self.INDEP_STRUC)
        self.brain_image_picker.Enable(mode == self.INDEP_STRUC and self.brain_image_picker.checkbox.IsChecked())

        # Only offer FSL_ANAT transform option if we are using FSL_ANAT
        sel = self.transform_ch.GetSelection()
        if mode == self.INDEP_STRUC: 
            if sel == self.TRANS_FSLANAT: sel = self.TRANS_MATRIX
            choices=2
        else: 
            choices=3
        self.transform_ch.Enable(False)
        self.transform_ch.Clear()
        self.transform_ch.AppendItems(self.transform_choices[:choices])
        self.transform_ch.SetSelection(sel)
        self.transform_ch.Enable(self.transform())

        self.transform_picker.Enable(self.transform() and self.transform_type() != self.TRANS_FSLANAT)

        TabPage.update(self)
        