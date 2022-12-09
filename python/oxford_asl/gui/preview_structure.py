"""
BASIL GUI for Oxford ASL - Preview of the structure of the data set

Copyright (C) 2020 University of Oxford
"""
import wx

class DataStructurePreview(wx.Panel):
    """
    Visual preview of the structure of an ASL data set
    """

    def __init__(self, parent, ntis, repeats, ntc, order, tagfirst):
        wx.Panel.__init__(self, parent, size=wx.Size(300, 150))
        self.SetBackgroundStyle(wx.BG_STYLE_CUSTOM)
        self.Bind(wx.EVT_SIZE, self._on_size)
        self.Bind(wx.EVT_PAINT, self._on_paint)
        self.ntis = ntis
        self.ntc = ntc
        self.repeats = repeats
        self.tagfirst = tagfirst
        self.order = order
        self.tis_name = "PLD"
        self.order_labels = {
            "r" : ("Repeat ", "R"),
            "t" : (self.tis_name, self.tis_name),
            "p" : {
                "tc" : (("Label", "Control"), ("L", "C")),
                "ct" : (("Control", "Label"), ("C", "L")),
            }
        }

        self.hfactor = 0.95
        self.vfactor = 0.95
        self.cols = {
            "r" : wx.Colour(128, 128, 255),
            "t" : wx.Colour(255, 128, 128),
            "p" : wx.Colour(128, 255, 128),
        }
        self.num = {"t" : ntis, "r" : repeats[0], "p" : ntc}

    def _on_size(self, event):
        event.Skip()
        self.Refresh()

    def _on_paint(self, _event):
        self.order_labels["t"] = (self.tis_name, self.tis_name)
        self.num = {"t" : self.ntis, "r" : self.repeats[0], "p" : self.ntc}

        if self.ntc == 1:
            order = self.order.replace("p", "")
        else:
            order = self.order

        width, height = self.GetClientSize()
        group_height = int(0.8*self.vfactor*height / len(order))
        group_width = int(self.hfactor*width)
        ox = int(width*(1-self.hfactor)/2)
        oy = int(height*(1-self.vfactor)/2)

        dc = wx.AutoBufferedPaintDC(self)
        dc.Clear()
        nvols = int(sum(self.repeats) * self.ntc)
        smallest_group_width = self._draw_groups(dc, order[::-1], ox, oy, group_width, group_height)
        self._centered_text(dc, "Input data volumes", ox+group_width/2, oy+0.9*height)
        self._centered_text(dc, "1", ox+smallest_group_width/2, oy+0.85*height)
        self._centered_text(dc, str(nvols), ox+group_width-smallest_group_width/2, oy+0.85*height)

    def _get_label(self, code, num, short):
        labels = self.order_labels[code]
        if isinstance(labels, dict):
            iaf = "tc" if self.tagfirst else "ct"
            labels = labels[iaf]
        label = labels[int(short)]
        if isinstance(label, tuple):
            return label[num]
        else:
            return label + str(int(num+1))

    def _centered_text(self, dc, text, x, y):
        text_size = dc.GetTextExtent(text)
        dc.DrawText(text, int(x-text_size.x/2), int(y-text_size.y/2))

    def _draw_groups(self, dc, groups, ox, oy, width, height, cont=False):
        smallest_width = width

        if groups:
            small = width < 150 # Heuristic
            group = groups[0]
            col = self.cols[group]
            if cont:
                # This 'group' is a continuation ellipsis
                rect = wx.Rect(ox, oy, width-1, height-1)
                dc.SetBrush(wx.Brush(col, wx.SOLID))
                dc.DrawRectangle(*rect.Get())
                text_size = dc.GetTextExtent("...")
                dc.DrawText("...", int(ox+width/2-text_size.x/2), int(oy+height/2-text_size.y/2))

                # Continuation ellipsis contains similar box for each group below it
                smallest_width = self._draw_groups(dc, groups[1:], ox, oy+height, width, height, cont=True)
            else:
                num = self.num[group]
                # Half the width of a normal box (full width of ellipsis box)
                box_width = int(width/min(2*num, 5))

                # Draw first
                label = self._get_label(group, 0, small)
                rect = wx.Rect(ox, oy, 2*box_width-1, height-1)
                dc.SetBrush(wx.Brush(col, wx.SOLID))
                dc.DrawRectangle(*rect.Get())
                text_size = dc.GetTextExtent(label)
                dc.DrawText(label, int(ox+box_width-text_size.x/2), int(oy+height/2-text_size.y/2))

                # Draw groups inside this group
                smallest_width = self._draw_groups(dc, groups[1:], ox, oy+height, 2*box_width, height)
                ox += 2*box_width

                # Draw ellipsis if required
                if num > 2:
                    smallest_width = self._draw_groups(dc, groups, ox, oy, box_width, height, cont=True)
                    ox += box_width

                # Draw last box if required
                if num > 1:
                    label = self._get_label(group, num-1, small)

                    rect = wx.Rect(ox, oy, 2*box_width-1, height-1)
                    dc.SetBrush(wx.Brush(col, wx.SOLID))
                    dc.DrawRectangle(*rect.Get())
                    text_size = dc.GetTextExtent(label)
                    dc.DrawText(label, int(ox+box_width-text_size.x/2), int(oy+height/2-text_size.y/2))

                    smallest_width = self._draw_groups(dc, groups[1:], ox, oy+height, 2*box_width, height)

        return smallest_width
