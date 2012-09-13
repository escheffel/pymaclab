'''
.. module:: queued
   :platform: Linux
   :synopsis: A collection of tools required for doing intelligent and dynamic DSGE model instance updating at runtime. The version in
              this module is for the queued updater behaviour.

.. moduleauthor:: Eric M. Scheffel <eric.scheffel@nottingham.edu.cn>


'''
from copy import deepcopy
from ..solvers.steadystate import ManualSteadyState

queue = []

class Updaters_Queued(object):
    def __init__(self):
        pass

class dicwrap_queued:
    def __init__(self,other,wrapobj_str,initlev):
        self.other = other
        self.queue = other.updaters_queued.queue
        self.wrapobj_str = wrapobj_str
        wrapobj = eval('other.'+wrapobj_str.split('.')[1])
        self.wrapobj = wrapobj
        self.initlev = initlev
        if wrapobj_str == 'self.vardic':
            if 'vardic' not in dir(other.updaters):
                self.wrapdic = deepcopy(other.vardic)
            else:
                self.wrapdic = other.updaters.vardic
        elif wrapobj_str == 'self.nlsubsdic':
            if 'nlsubsdic' not in dir(other.updaters):
                self.wrapdic = other.nlsubsdic.copy()
            else:
                self.wrapdic = other.updaters.nlsubsdic
        elif wrapobj_str == 'self.paramdic':
            if 'paramdic' not in dir(other.updaters):
                self.wrapdic = other.paramdic.copy()
            else:
                self.wrapdic = other.updaters.paramdic
            
    def __getattr__(self,attrname):
        return getattr(self.wrapdic,attrname)

    def __setitem__(self,key,value):
        other = self.other
        initlev = self.initlev
        wrapobj_str = self.wrapobj_str
        wrapobj = self.wrapobj
        # ...and assign before test of inequality
        wrapobj[key] = value
        if self.wrapdic != wrapobj:
            self.wrapdic.update(wrapobj)
            if wrapobj_str == 'self.nlsubsdic':
                for i1,elem in enumerate(other.nlsubs_raw1):
                    other.nlsubs_raw1[i1][1] = deepcopy(self.nlsubsdic[other.nlsubs_raw1[i1][0]])
                if 'self.nlsubsdic' not in self.queue: self.queue.append('self.nlsubsdic')
            elif wrapobj_str == 'self.paramdic':
                if 'self.paramdic' not in self.queue: self.queue.append('self.paramdic')
            elif wrapobj_str == 'self.vardic':
                if 'self.vardic' not in self.queue: self.queue.append('self.vardic')           

    def __getitem__(self,key):
        return self.wrapdic[key]

    def update(self,dico):
        self.updic = dico
        other = self.other
        initlev = self.initlev
        wrapobj_str = self.wrapobj_str
        wrapobj = self.wrapobj
        # ...and update befor test of inequality
        wrapobj.update(dico)
        # Check if new keys are already present in wrapdic
        for keyo in dico.keys():
            if keyo not in self.wrapdic.keys():
                print "ERROR: You can only update existing keys, not introduce new ones."
                return
        # Check if any key's value has been updated relative to wrapdic
        if self.wrapdic != wrapobj:
            self.wrapdic.update(wrapobj)
            if wrapobj_str == 'self.nlsubsdic':
                for i1,elem in enumerate(other.nlsubs_raw1):
                    other.nlsubs_raw1[i1][1] = deepcopy(self.nlsubsdic[other.nlsubs_raw1[i1][0]])
                if 'self.nlsubsdic' not in self.queue: self.queue.append('self.nlsubsdic')
            elif wrapobj_str == 'self.paramdic':
                if 'self.paramdic' not in self.queue: self.queue.append('self.paramdic')
            elif wrapobj_str == 'self.vardic':
                if 'self.vardic' not in self.queue: self.queue.append('self.vardic')           


    def __repr__(self):
        return self.wrapdic.__repr__()
    def __str__(self):
        return self.wrapdic.__str__()

class listwrapk:
    def __init__(self,other,wrapobj,initlev):
        self.other = other
        self.wrapobj = wrapobj
        self.initlev = initlev
        self.wrapli = deepcopy(wrapobj)
        self.queue = other.updaters_queued.queue
            
    def __getattr__(self,attrname):
        return getattr(self.wrapli,attrname)

    def __setslice__(self,ind1,ind2,into):
        other = self.other
        wrapobj = self.wrapobj
        initlev = self.initlev
        lengo = len(self.wrapli)
        if ind2 >= lengo:
            print "ERROR: Assignment out of bounds of original list"
            return
        ##### THE INITS #####################
        #other.init1()
        if self.wrapli[ind1:ind2] != into:
            self.wrapli[ind1:ind2] = into
            self.queue.append('self.vardic')
    
    def __setitem__(self,ind,into):
        other = self.other
        initlev = self.initlev
        wrapobj = self.wrapobj
        lengo = len(self.wrapli)
        if ind >= lengo:
            print "ERROR: Assignment out of bounds of original list"
            return
        ##### THE INITS #####################
        #other.init1()
        if self.wrapli[ind] != into:
            self.wrapli[ind] = into
            if 'self.vardic' not in self.queue: self.queue.append('self.vardic')
            
                  

    def __getitem__(self,ind):
        lengo = len(self.wrapli)
        if ind >= lengo:
            print "ERROR: Assignment out of bounds of original list"
            return
        else:
            return self.wrapli[ind]

    def __repr__(self):
        return self.wrapli.__repr__()
    def __str__(self):
        return self.wrapli.__str__()

class dicwrapk:
    def __init__(self,other,wrapobj,initlev):
        self.other = other
        self.wrapobj = wrapobj
        self.initlev = initlev
        self.wrapdic = deepcopy(wrapobj)
        self.queue = other.updaters_queued.queue

    def __getattr__(self,attrname):
        return getattr(self.wrapdic,attrname)

    def __setitem__(self,key,value):
        other = self.other
        initlev = self.initlev
        wrapobj = self.wrapobj
        wrapobj[key] = value
        # Test if the dictionary has changed relative to self.wrapdic
        if self.wrapdic != wrapobj:
            self.wrapdic[key] = deepcopy(wrapobj[key])
            ##### THE INITS #####################
            #other.init1()
            if 'self.vardic' not in self.queue: self.queue.append('self.vardic')


    def __getitem__(self,key):
        return self.wrapdic[key]

    def update(self,dico):
        self.updic = dico
        other = self.other
        initlev = self.initlev
        wrapobj = self.wrapobj
        wrapobj.update(dico)
        # Check if new keys are already present in wrapdic
        for keyo in dico.keys():
            if keyo not in self.wrapdic.keys():
                print "ERROR: You can only update existing keys, not introduce new ones."
                return
        # Check if any key's value has been updated relative to wrapdic
        if self.wrapdic != wrapobj:
            self.wrapdic.update(wrapobj)
            ##### THE INITS #####################
            #other.init1()
            if 'self.vardic' not in self.queue: self.queue.append('self.vardic')


    def __repr__(self):
        return self.wrapdic.__repr__()
    def __str__(self):
        return self.wrapdic.__str__()    

class dicwrap_deep_queued:
    def __init__(self,other,wrapobj_str,initlev):
        self.other = other
        self.wrapobj_str = wrapobj_str
        self.initlev = initlev
        self.queue = other.updaters_queued.queue
        # Create the wrapobj using the passed string
        self.wrapobj = eval('other.'+wrapobj_str.split('.')[1])
        if wrapobj_str == 'self.vardic':
            self.wrapdic = deepcopy(self.wrapobj)
            for keyo in self.wrapdic.keys():
                self.wrapdic[keyo] = dicwrapk(other,self.wrapdic[keyo],initlev)
                for keyo2 in self.wrapdic[keyo].keys():
                    self.wrapdic[keyo][keyo2] = dicwrapk(other,self.wrapdic[keyo][keyo2],initlev)
                    for i1,elem in enumerate(self.wrapdic[keyo][keyo2]):
                        self.wrapdic[keyo][keyo2][i1] = listwrapk(other,self.wrapdic[keyo][keyo2][i1],initlev)

    def __getattr__(self,attrname):
        return getattr(self.wrapdic,attrname)
                       
    def __setitem__(self,key,value):
        other = self.other
        initlev = self.initlev
        wrapobj_str = self.wrapobj_str
        wrapobj = self.wrapobj
        wrapobj[key] = value
        # Test if the dictionary has changed relative to self.wrapdic
        if self.wrapdic != wrapobj:
            self.wrapdic[key] = value
            ##### THE INITS #####################
            #other.init1()
            if wrapobj_str == 'self.vardic':
                if 'self.vardic' not in self.queue: self.queue.append('self.vardic')

            other.init1a()
            if wrapobj_str == 'self.nlsubsdic':
                # not a deep dic
                pass

            other.init1b()
            if wrapobj_str == 'self.paramdic':
                # not a deep dic
                pass


    def __update__(self,dico):
        other = self.other
        initlev = self.initlev
        wrapobj_str = self.wrapobj_str
        wrapobj = self.wrapobj
        wrapobj.update(dico)
        # Test if the dictionary has changed relative to self.wrapdic
        if self.wrapdic != wrapobj:
            self.wrapdic.update(dico)
            ##### THE INITS #####################
            #other.init1()
            if wrapobj_str == 'self.vardic':
                if 'self.vardic' not in self.queue: self.queue.append('self.vardic')



class listwrap_queued:
    def __init__(self,other,wrapobj_str,initlev):
        self.other = other
        self.queue = other.updaters_queued.queue
        self.wrapobj_str = wrapobj_str
        wrapobj = eval('other.'+wrapobj_str.split('.')[1])
        self.wrapobj = wrapobj
        self.initlev = initlev
        self.wrapli = deepcopy(wrapobj)

    def __setslice__(self,ind1,ind2,into):
        other = self.other
        wrapobj_str = self.wrapobj_str
        wrapobj = self.wrapobj
        initlev = self.initlev
        lengo = len(self.wrapli)
        if ind2 >= lengo:
            print "ERROR: Assignment out of bounds of original list"
            return

        if self.wrapli[ind1:ind2] != into and wrapobj_str == 'self.foceqs':
            self.wrapli[ind1:ind2] = into
            if 'self.foceqs' not in self.queue: self.queue.append('self.foceqs')
        elif self.wrapli[ind1:ind2] != into and wrapobj_str == 'self.manss_sys':
            self.wrapli[ind1:ind2] = into
            if 'self.manss_sys' not in self.queue: self.queue.append('self.manss_sys')
        elif self.wrapli[ind1:ind2] != into and wrapobj_str == 'self.ssys_list':
            self.wrapli[ind1:ind2] = into
            if 'self.ssys_list' not in self.queue: self.queue.append('self.ssys_list')
    
    def __setitem__(self,ind,into):
        other = self.other
        wrapobj_str = self.wrapobj_str
        wrapobj = self.wrapobj
        initlev = self.initlev
        lengo = len(self.wrapli)
        if ind >= lengo:
            print "ERROR: Assignment out of bounds of original list"
            return

        if self.wrapli[ind] != into and wrapobj_str == 'self.foceqs':
            self.wrapli[ind] = into
            if 'self.foceqs' not in self.queue: self.queue.append('self.foceqs')
        elif self.wrapli[ind] != into and wrapobj_str == 'self.manss_sys':
            self.wrapli[ind] = into
            if 'self.manss_sys' not in self.queue: self.queue.append('self.manss_sys')
        elif self.wrapli[ind] != into and wrapobj_str == 'self.ssys_list':
            self.wrapli[ind] = into
            if 'self.ssys_list' not in self.queue: self.queue.append('self.ssys_list')

    def __getitem__(self,ind):
        lengo = len(self.wrapli)
        if ind >= lengo:
            print "ERROR: Assignment out of bounds of original list"
            return
        else:
            return self.wrapli[ind]

    def __repr__(self):
        return self.wrapli.__repr__()
    def __str__(self):
        return self.wrapli.__str__()


class matwrap_queued:
    def __init__(self,other,wrapobj_str,initlev):
        self.other = other
        self.queue = other.updaters_queued.queue
        self.wrapobj_str = wrapobj_str
        self.initlev = initlev
        if wrapobj_str == 'self.sigma':
            self.wrapmat = deepcopy(other.sigma)
            
    def __getattr__(self,attrname):
        return getattr(self.wrapmat,attrname)
    
    def __setitem__(self,ind,into):
        # ind is a tuple and into is (should be) a float
        other = self.other
        wrapob_str = self.wrapobj_str
        initlev = self.initlev
        if type(into) != type(1.111): into = float(into)
        if self.wrapmat[ind[0],ind[1]] != into and wrapob_str == 'self.sigma':
            self.wrapmat[ind[0],ind[1]] = into
            other.sigma[ind[0],ind[1]] = into
            if 'self.sigma' not in self.queue: self.queue.append('self.sigma')

   
class Process_Queue(object):
    def __init__(self,other=None):
        self.other = other
        self.queue = other.updaters_queued.queue
        self.initlev = other._initlev
        # The vardic
        self.vardic = other.updaters_queued.vardic
        # The nlsubsdic
        self.nlsubsdic = other.updaters_queued.nlsubsdic
        # The paramdic
        self.paramdic = other.updaters_queued.paramdic
        # The foceqs
        self.foceqs = other.updaters_queued.foceqs
        # The manss_sys
        if 'manss_sys' in dir(other.updaters_queued):
            self.manss_sys = other.updaters_queued.manss_sys
        # The ssys_list
        if 'ssys_list' in dir(other.updaters_queued):
            self.ssys_list = other.updaters_queued.ssys_list
        # The sigma
        self.sigma = other.updaters_queued.sigma
        
    def __call__(self):
        queue = self.queue
        other = self.other
        initlev = self.initlev
        mesg = other._mesg
        if mesg:
            print "You have UPDATED the following items which are now PROCESSED:"
            for elem in queue:
                print elem
            print '================================================================'
        ##### THE INITS #####################
        other.init1()
        if 'self.vardic' in queue:
            other.vardic = deepcopy(self.vardic)

        other.init1a()
        if 'self.nlsubsdic' in queue:
            for i1,elem in enumerate(other.nlsubs_raw1):
                other.nlsubs_raw1[i1][1] = deepcopy(self.nlsubsdic[other.nlsubs_raw1[i1][0]])
            other.nlsubsdic = self.nlsubsdic.copy()

        other.init1b()
        if 'self.paramdic' in queue:
            other.paramdic = self.paramdic.copy()
        
        other.init1c()
        if 'self.foceqs' in queue:
            other.foceqs = deepcopy(self.foceqs)
        if 'self.manss_sys' in queue:
            other.manss_sys = deepcopy(self.manss_sys)
        if 'self.ssys_list' in queue:
            other.ssys_list = deepcopy(self.ssys_list)

        # Prepare DSGE model instance for manual SS computation
        other.init2()
        if initlev == 0:
            other.init_out()

        # Solve for SS automatically
        other.init3()
        if initlev == 1:
            other.init_out() 

        other.init4()
        if 'self.sigma' in queue:
            other.sigma = deepcopy(self.sigma)

        # Solve model dynamically            
        other.init5()
        if initlev == 2:
            other.init_out()
            
        return
