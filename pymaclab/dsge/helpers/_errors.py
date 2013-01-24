import helpers as HLP

# ERROR CLASSES
class Errbase(Exception):
	pass

class Uhlerr(Errbase):
	pass

class Dataerr(Errbase):
	pass

class TS_err(Errbase):
	def __init__(self,inmes):
		self.inmes = inmes
		now = HLP.now_is()
		if self.inmes == 'Singular Matrix':
			self.outmes = 'Singular Matrix XX!! Cannot solve for beta'
			print str(now) + '   ' + self.outmes
			
# ERROR CLASSES OUT