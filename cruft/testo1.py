import wx
import wx.gizmos as gizmos

class AlarmsWindow(wx.MDIChildFrame):
	def __init__(self, parent, title, size, pos):
		wx.MDIChildFrame.__init__(self, parent, -1, title, size = size, 
pos =pos)
		self.tree = gizmos.TreeListCtrl(self, -1, style = 
wx.TR_DEFAULT_STYLE| wx.TR_FULL_ROW_HIGHLIGHT)
	# create some columns
		self.tree.AddColumn("Connection")
		self.tree.AddColumn("Alarm")
		self.tree.AddColumn("Value")
		self.tree.SetMainColumn(0)

		self.root = self.tree.AddRoot("Connections")
		self.tree.Expand(self.root)
		self.tree.GetMainWindow().Bind(wx.EVT_RIGHT_UP, self.OnRightUp)
		self.Show()

		child = self.tree.AppendItem(self.root, 'name')
		self.tree.SetItemText(child, 'name')
		child2= self.tree.AppendItem(child,'name2')
##        x = self.tree.AppendItem(child2, "")
		x = self.tree.AppendItem(child2, "XXX")
		self.tree.SetItemText(x, 'alarm', 1)
		self.tree.SetItemText(x, 'value', 2)

	def OnRightUp(self, evt):
		pass

class MDIFrame(wx.MDIParentFrame):
	def __init__(self):
		wx.MDIParentFrame.__init__(self, None, -1, "MDI Parent", size 
=(600, 400))
		child = AlarmsWindow(self, "Alarm", (400, 300), (0, 0))


if __name__=='__main__':
	app = wx.PySimpleApp()
	frame = MDIFrame()
	frame.Show()
	app.MainLoop()
