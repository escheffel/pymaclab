"""
COLLECTION OF SUPPORTING FUNCTIONS AND CLASSES
"""
from copy import deepcopy
from ..solvers.steadystate import ManualSteadyState

class Updaters(object):
    def __init__(self):
        pass

class dicwrap:
    def __init__(self,other,wrapobj_str,initlev):
        self.other = other
	self.wrapobj_str = wrapobj_str
        self.initlev = initlev
	if wrapobj_str == 'self.nlsubsdic':
	    self.wrapdic = other.nlsubsdic
	elif wrapobj_str == 'self.paramdic':
	    self.wrapdic = other.paramdic
    def __getattr__(self,attrname):
        return getattr(self.wrapdic,attrname)

    def __setitem__(self,key,value):
        other = self.other
        initlev = self.initlev
	wrapobj_str = self.wrapobj_str
        if self.wrapdic[key] != value:
            self.wrapdic[key] = value
	    ##### THE INITS #####################
	    other.init1()
	    if wrapobj_str == 'self.nlsubsdic':
		other.nlsubsdic = deepcopy(self.wrapdic)
		for i1,elem in enumerate(other.nlsubs_raw1):
		    other.nlsubs_raw1[i1][1] = self.wrapdic[other.nlsubs_raw1[i1][0]]
		other.nlsubsdic = deepcopy(self.wrapdic)

	    other.init1a()
	    if wrapobj_str == 'self.paramdic':
		other.paramdic = deepcopy(self.wrapdic)

	    other.init1b()

	    # Prepare DSGE model instance for manual SS computation
	    other.init2()
	    if initlev == 0:
		    other.init_out()

	    # Solve for SS automatically
	    other.init3()
	    if initlev == 1:
		    other.init_out()

	    # Solve model dynamically
	    other.init4()
	    if initlev == 2:
		    other.init_out()

    def __getitem__(self,key):
        return self.wrapdic[key]

    def update(self,dico):
	self.updic = dico
	other = self.other
	initlev = self.initlev
	wrapobj_str = self.wrapobj_str
	# Check if new keys are already present in wrapdic
	for keyo in dico.keys():
	    if keyo not in self.wrapdic.keys():
		print "ERROR: You can only update existing keys, not introduce new ones."
		return
	# Check if any key's value has been updated relative to wrapdic
	if any([True if dico[keyo]!= self.wrapdic[keyo] else False for keyo in dico]):
	    self.wrapdic.update(dico)
	    ##### THE INITS #####################
	    other.init1()
	    if wrapobj_str == 'self.nlsubsdic':
		other.nlsubsdic = deepcopy(self.wrapdic)
		for i1,elem in enumerate(other.nlsubs_raw1):
		    other.nlsubs_raw1[i1][1] = self.wrapdic[other.nlsubs_raw1[i1][0]]
		other.nlsubsdic = deepcopy(self.wrapdic)

	    other.init1a()
	    if wrapobj_str == 'self.paramdic':
		other.paramdic = deepcopy(self.wrapdic)

	    other.init1b()

	    # Prepare DSGE model instance for manual SS computation
	    other.init2()
	    if initlev == 0:
		    other.init_out()

	    # Solve for SS automatically
	    other.init3()
	    if initlev == 1:
		    other.init_out()

	    # Solve model dynamically
	    other.init4()
	    if initlev == 2:
		    other.init_out()
	    
    
    def __repr__(self):
        return self.wrapdic.__repr__()
    def __str__(self):
        return self.wrapdic.__str__()
