import wx

class MyFrame(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, wx.DefaultPosition)

        sizer = wx.GridBagSizer(9, 9)
        sizer.Add(wx.Button(self,-1, "Button"), (0, 0), wx.DefaultSpan,  wx.ALL, 5)
        sizer.Add(wx.Button(self,-1, "Button"), (1, 1), (1,7), wx.EXPAND)
        sizer.Add(wx.Button(self,-1, "Button"), (6, 6), (3,3), wx.EXPAND)
        sizer.Add(wx.Button(self,-1, "Button"), (3, 0), (1,1), wx.ALIGN_CENTER)
        sizer.Add(wx.Button(self,-1, "Button"), (4, 0), (1,1), wx.ALIGN_LEFT)
        sizer.Add(wx.Button(self,-1, "Button"), (5, 0), (1,1), wx.ALIGN_RIGHT)
        sizer.AddGrowableRow(6)
        sizer.AddGrowableCol(6)

        self.SetSizerAndFit(sizer)
        self.Centre()

class MyApp(wx.App):
    def OnInit(self):
        frame = MyFrame(None, -1, "wxgridbagsizer.py")
        frame.Show(True)
        self.SetTopWindow(frame)
        return True

app = MyApp(0)
app.MainLoop()