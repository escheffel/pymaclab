'''
.. module:: one_off
   :platform: Linux
   :synopsis: A collection of tools required for doing intelligent and dynamic DSGE model instance updating at runtime. The version in
              this module is for the one-off updater behaviour.

.. moduleauthor:: Eric M. Scheffel <eric.scheffel@nottingham.edu.cn>


'''
from copy import deepcopy
from ..solvers.steadystate import ManualSteadyState

class Updaters(object):
    def __init__(self):
        pass

class dicwrap:
    def __init__(self,other,wrapobj_str,initlev):
        self.other = other
        self.wrapobj_str = wrapobj_str
        wrapobj = eval('other.'+wrapobj_str.split('.')[1])
        self.wrapobj = wrapobj
        self.initlev = initlev
        if wrapobj_str == 'self.vardic':
            if 'vardic' not in dir(other.updaters_queued):
                self.wrapdic = deepcopy(other.vardic)
            else:
                self.wrapdic = other.updaters_queued.vardic
        elif wrapobj_str == 'self.nlsubsdic':
            if 'nlsubsdic' not in dir(other.updaters_queued):
                self.wrapdic = other.nlsubsdic.copy()
            else:
                self.wrapdic = other.updaters_queued.nlsubsdic
        elif wrapobj_str == 'self.paramdic':
            if 'paramdic' not in dir(other.updaters_queued):
                self.wrapdic = other.paramdic.copy()
            else:
                self.wrapdic = other.updaters_queued.paramdic

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
            queue1 = deepcopy(other.updaters_queued.queue)
            self.wrapdic.update(wrapobj)
            queue2 = deepcopy(other.updaters_queued.queue)            
            ##### THE INITS #####################
            other.init1()
            # Remove from queue
            if queue1 != queue2:
                for elem in queue2:
                    if elem not in queue1: other.updaters_queued.queue.remove(elem)
            if wrapobj_str == 'self.vardic':
                wrapobj.update(self.wrapdic)

            other.init1a()
            if wrapobj_str == 'self.nlsubsdic':
                for i1,elem in enumerate(other.nlsubs_raw1):
                    other.nlsubs_raw1[i1][1] = self.wrapdic[other.nlsubs_raw1[i1][0]]
                other.nlsubsdic.update(self.wrapdic)

            other.init1b()
            if wrapobj_str == 'self.paramdic':
                other.paramdic.update(self.wrapdic)
            
            other.init1c()

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
            other.init5()
            if initlev == 2:
                other.init_out()

    def __getitem__(self,key):
        return self.wrapdic[key]

    def update(self,dico):
        self.updic = dico
        other = self.other
        initlev = self.initlev
        wrapobj_str = self.wrapobj_str
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
            other.init1()
            if wrapobj_str == 'self.vardic':
                other.vardic = deepcopy(self.wrapdic)

            other.init1a()
            if wrapobj_str == 'self.nlsubsdic':
                for i1,elem in enumerate(other.nlsubs_raw1):
                    other.nlsubs_raw1[i1][1] = self.wrapdic[other.nlsubs_raw1[i1][0]]
                other.nlsubsdic.update(self.wrapdic)

            other.init1b()
            if wrapobj_str == 'self.paramdic':
                other.paramdic.update(self.wrapdic)
            
            other.init1c()

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
            other.init5()
            if initlev == 2:
                other.init_out()


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
            wrapobj[ind1:ind2] = into
            
            other.init1a()
        
            other.init1b()
        
            other.init1c()

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
            other.init5()
            if initlev == 2:
                other.init_out() 
    
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
            wrapobj[ind] = into
            
            other.init1a()
            
            other.init1b()
            
            other.init1c()            
            
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
            other.init5()
            if initlev == 2:
                other.init_out()             
                  

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

    def __getattr__(self,attrname):
        return getattr(self.wrapdic,attrname)

    def __setitem__(self,key,value):
        other = self.other
        initlev = self.initlev
        wrapobj = self.wrapobj
        wrapobj[key] = value
        # Test if the dictionary has changed relative to self.wrapdic
        if self.wrapdic != wrapobj:
            self.wrapdic[key] = value
            ##### THE INITS #####################
            #other.init1()
            other.vardic.upate(self.wrapdic)

            other.init1a()

            other.init1b()
            
            other.init1c()

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
            other.init5()
            if initlev == 2:
                other.init_out()

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
            self.wrapdic.update(dico)
            ##### THE INITS #####################
            #other.init1()
            other.vardic = deepcopy(self.wrapdic)

            other.init1a()

            other.init1b()
            
            other.init1c()

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
            other.init5()
            if initlev == 2:
                other.init_out()


    def __repr__(self):
        return self.wrapdic.__repr__()
    def __str__(self):
        return self.wrapdic.__str__()
   
class dicwrap_deep:
    def __init__(self,other,wrapobj_str,initlev):
        self.other = other
        self.wrapobj_str = wrapobj_str
        self.initlev = initlev
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
        wrapobj = self.wrapobj
        wrapobj[key] = value
        # Test if the dictionary has changed relative to self.wrapdic
        if self.wrapdic != wrapobj:
            self.wrapdic[key] = value
            ##### THE INITS #####################
            #other.init1()
            if wrapobj_str == 'self.vardic':
                other.paramdic.upate(self.wrapdic)

            other.init1a()
            if wrapobj_str == 'self.nlsubsdic':
                # not a deep dic
                pass

            other.init1b()
            if wrapobj_str == 'self.paramdic':
                # not a deep dic
                pass
            
            other.init1c()

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
            other.init5()
            if initlev == 2:
                other.init_out()


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
                other.paramdic.upate(self.wrapdic)

            other.init1a()

            other.init1b()
            
            other.init1c()

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
            other.init5()
            if initlev == 2:
                other.init_out()
                


class listwrap:
    def __init__(self,other,wrapobj_str,initlev):
        self.other = other
        self.wrapobj_str = wrapobj_str
        wrapobj = eval('other.'+wrapobj_str.split('.')[1])
        self.wrapobj = wrapobj
        self.initlev = initlev
        self.wrapli = deepcopy(wrapobj)

            
    def __getattr__(self,attrname):
        return getattr(self.wrapli,attrname)

    def __setslice__(self,ind1,ind2,into):
        other = self.other
        wrapobj_str = self.wrapobj_str
        wrapobj = eval('other.'+wrapobj_str.split('.')[1])
        initlev = self.initlev
        do_rest = False
        lengo = len(self.wrapli)
        if ind2 >= lengo:
            print "ERROR: Assignment out of bounds of original list"
            return
        ##### THE INITS #####################
        #other.init1()
        #other.init1a()
        #other.init1b()
        #other.init1c()
        if self.wrapli[ind1:ind2] != into:
            self.wrapli[ind1:ind2] = into
            wrapobj[ind1:ind2] = into
            do_rest = True

        if do_rest:
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
            other.init5()
            if initlev == 2:
                other.init_out() 
    
    def __setitem__(self,ind,into):
        other = self.other
        wrapobj_str = self.wrapobj_str
        wrapobj = eval('other.'+wrapobj_str.split('.')[1])
        initlev = self.initlev
        do_rest = False
        lengo = len(self.wrapli)
        if ind >= lengo:
            print "ERROR: Assignment out of bounds of original list"
            return
        ##### THE INITS #####################
        #other.init1()
        #other.init1a()
        #other.init1b()
        #other.init1c()
        if self.wrapli[ind] != into:
            self.wrapli[ind] = into
            wrapobj[ind] = into
            do_rest = True

        if do_rest:
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
            other.init5()
            if initlev == 2:
                other.init_out()           

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
    

class matwrap:
    def __init__(self,other,wrapobj_str,initlev):
        self.other = other
        self.wrapobj_str = wrapobj_str
        self.initlev = initlev
        if wrapobj_str == 'self.sigma':
            self.wrapmat = other.sigma
            
    def __getattr__(self,attrname):
        return getattr(self.wrapmat,attrname)
    
    def __setitem__(self,ind,into):
        other = self.other
        wrapob_str = self.wrapobj_str
        initlev = self.initlev
        ##### THE INITS #####################
        other.init1()
        other.init1a()
        other.init1b()
        other.init1c()

        # Prepare DSGE model instance for manual SS computation
        other.init2()
        if initlev == 0:
            other.init_out()

        # Solve for SS automatically
        other.init3()
        if initlev == 1:
            other.init_out()

        other.init4()
        if self.wrapmat[ind[0],ind[1]] != into and wrapob_str == 'self.sigma':
            self.wrapmat[ind[0],ind[1]] = into
            other.sigma[ind[0],ind[1]] = into
            
        # Solve model dynamically
        other.init5()
        if initlev == 2:
            other.init_out()  