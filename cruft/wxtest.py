import wx
import os
ID_ABOUT=101
ID_OPEN=102
ID_BUTTON1=110
ID_EXIT=200
class MainWindow(wx.Frame):
    def __init__(self,parent,id,title):
        self.dirname=''
        wx.Frame.__init__(self,parent,wx.ID_ANY, title)
        self.control = wx.TextCtrl(self, 1, style=wx.TE_MULTILINE)
        self.CreateStatusBar() # A Statusbar in the bottom of the window
        # Setting up the menu.
        filemenu= wx.Menu()
        filemenu.Append(ID_OPEN, "&Open"," Open a file to edit")
        filemenu.AppendSeparator()
        filemenu.Append(ID_ABOUT, "&About"," Information about this program")
        filemenu.AppendSeparator()
        filemenu.Append(ID_EXIT,"E&xit"," Terminate the program")
        # Creating the menubar.
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu,"&File") # Adding the "filemenu" to the MenuBar
        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
        wx.EVT_MENU(self, ID_ABOUT, self.OnAbout)
        wx.EVT_MENU(self, ID_EXIT, self.OnExit)
        wx.EVT_MENU(self, ID_OPEN, self.OnOpen)
        self.sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        self.buttons=[]
        for i in range(0,6):
            self.buttons.append(wx.Button(self, ID_BUTTON1+i, "Button &"+`i`))
            self.sizer2.Add(self.buttons[i],1,wx.EXPAND)
        # Use some sizers to see layout options
        self.sizer=wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.control,1,wx.EXPAND)
        self.sizer.Add(self.sizer2,0,wx.EXPAND)
        #Layout sizers
        self.SetSizer(self.sizer)
        self.SetAutoLayout(1)
        self.sizer.Fit(self)
        self.Show(1)
    def OnAbout(self,e):
        d= wx.MessageDialog( self, " A sample editor \n"
                            " in wxPython","About Sample Editor", wx.OK)
                            # Create a message dialog box
        d.ShowModal() # Shows it
        d.Destroy() # finally destroy it when finished.
    def OnExit(self,e):
        self.Close(True)  # Close the frame.
    def OnOpen(self,e):
        """ Open a file"""
        dlg = wx.FileDialog(self, "Choose a file", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename=dlg.GetFilename()
            self.dirname=dlg.GetDirectory()
            f=open(os.path.join(self.dirname, self.filename),'r')
            self.control.SetValue(f.read())
            f.close()
        dlg.Destroy()
app = wx.PySimpleApp()
frame = MainWindow(None, -1, "Sample editor")
app.MainLoop()