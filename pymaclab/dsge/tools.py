"""
COLLECTION OF SUPPORTING FUNCTIONS AND CLASSES
"""
from solvers.steadystate import ManualSteadyState

class dicwrap:
    def __init__(self,other,initlev):
        self.other = other
        self.initlev = initlev
        self.wrapdic = other.paramdic
    def __getattr__(self,attrname):
        return getattr(self.wrapdic,attrname)

    def __setitem__(self,key,value):
        other = self.other
        initlev = self.initlev
        nreg = self.nreg
        if self.wrapdic[key] != value:
            self.wrapdic[key] = value
            ################## STEADY STATE CALCULATIONS !!! ##############
            # Delete the models previous results
            del other.modsolvers
            del other.sssolvers
            other.sssolvers = SSsolvers()
            # Reset switch back to zero
            other.switches['ss_suc'] = ['0','0']
            # Solve for steady-state using fsolve
#            if sum([nreg.search(x)!=None for x in other.txtpars.secs['ssm'][0]]) == 0:
            if any([False if 'None' in x else True for x in secs['manualss'][0]]):
                intup = (other.ssys_list,other.ssidic,other.paramdic)
                other.sssolvers.fsolve = Fsolve(intup)
                other.sssolvers.fsolve.solve()
                if other.sssolvers.fsolve.ier == 1:
                    nsstate = other.sssolvers.fsolve.fsout
                    mreg = RE.compile('(.*?)_bar')
                    mreg2 = RE.compile('(.*?)_F|B.*?_bar')
                    for key in nsstate.keys():
                        ma = mreg.search(key)
                        if ma:
                            vari = ma.group(1)
                            for okey in other.sstate:
                                ma2 = mreg2.search(okey)
                                if ma2:
                                    vari2 = ma2.group(1)
                                    if vari == vari2:
                                        other.sstate[okey] = nsstate[key]
                    other.sstate.update(nsstate)
                    other.numssdic.update(nsstate)
                    other.switches['ss_suc'] = ['1','1']
                else:
                    other.switches['ss_suc'] = ['1','0']
            # Solve for steady-state using manss
#            if sum([nreg.search(x)!=None for x in other.txtpars.secs['sss'][0]]) == 0:
            if any([False if 'None' in x else True for x in secs['closedformss'][0]]):
                if other.switches['ss_suc'] == ['1','1']:
                    alldic = {}
                    alldic.update(other.sstate)
                    alldic.update(other.paramdic)
                    intup = (other.manss_sys,alldic)
                    other.sssolvers.manss = ManualSteadyState(intup)
                    other.sssolvers.manss.solve()
                    nsstate = other.sssolvers.manss.sstate
                    mreg = RE.compile('(.*?)_bar')
                    mreg2 = RE.compile('(.*?)_F|B.*?_bar')
                    for key in nsstate.keys():
                        ma = mreg.search(key)
                        if ma:
                            vari = ma.group(1)
                            for okey in other.sstate:
                                ma2 = mreg2.search(okey)
                                if ma2:
                                    vari2 = ma2.group(1)
                                    if vari == vari2:
                                        other.sstate[okey] = nsstate[key]
                    other.sstate.update(nsstate)
                else:
                    intup = (other.manss_sys,other.paramdic)
                    other.sssolvers.manss = ManualSteadyState(intup)
                    other.sssolvers.manss.solve()
                    nsstate = other.sssolvers.manss.sstate
                    mreg = RE.compile('(.*?)_bar')
                    mreg2 = RE.compile('(.*?)_F|B.*?_bar')
                    for key in nsstate.keys():
                        ma = mreg.search(key)
                        if ma:
                            vari = ma.group(1)
                            for okey in other.sstate:
                                ma2 = mreg2.search(okey)
                                if ma2:
                                    vari2 = ma2.group(1)
                                    if vari == vari2:
                                        other.sstate[okey] = nsstate[key]
                    other.sstate.update(nsstate)
            if initlev == 1: return
    
            # Open the model solution tree branch
            other.modsolvers = MODsolvers()
######################## LINEAR METHODS !!! ############################
#            if sum([nreg.search(x)!=None for x in other.txtpars.secs['modeq'][0]]) == 0:
            if any([False if 'None' in x else True for x in secs['modeq'][0]]):
                # Open the matlab Uhlig object
                intup = ((other.nendo,other.ncon,other.nexo),
                     other.eqindx,
                     other.vreg,
                     other.llsys_list,
                     other.diffli1,
                     other.diffli2,
                     sess1,
                     other.vardic)
                other.modsolvers.matuhlig = MatUhlig(intup)
                # Open the native Uhlig object
                intup = ((other.nendo,other.ncon,other.nexo),
                     other.eqindx,
                     other.vreg,
                     other.llsys_list,
                     other.diffli1,
                     other.diffli2,
                     sess1)
                other.modsolvers.pyuhlig = PyUhlig(intup)
                # Open the matlab Klein object
                intup = ((other.nendo,other.ncon,other.nexo),
                     other.eqindx,
                     other.vreg,
                     other.llsys_list,
                     other.diffli1,
                     other.diffli2,
                     sess1)
                other.modsolvers.matklein = MatKlein(intup)
                # Open the Fortran Klein object
                intup = ((other.nendo,other.ncon,other.nexo),
                     other.eqindx,
                     other.vreg,
                     other.llsys_list,
                     other.diffli1,
                     other.diffli2,
                     sess1)
                other.modsolvers.forklein = ForKlein(intup)
################## 1ST-ORDER NON-LINEAR METHODS !!! ##################
#            if sum([nreg.search(x)!=None for x in other.txtpars.secs['focs'][0]]) == 0:
            if any([False if 'None' in x else True for x in secs['focs'][0]]):
    
                # First, create the Jacobian and (possibly-->mk_hessian==True?) Hessian
                if use_anaderiv:
                    if ncpus > 1 and mk_hessian:
                        other.mkjaheppn()
                    elif ncpus > 1 and not mk_hessian:
                        other.mkjaheppn()
                    else:
                        other.mkjahen()
                else:
                    other.mkjahenmat()
    
                # Open the MatWood object
                intup = (other.jAA,other.jBB,
                     other.nexo,other.ncon,
                     other.nendo,sess1)
                other.modsolvers.matwood = MatWood(intup)
                # Open the Fortran KleinD object
                if 'nlsubsys' in dir(other):
                    intup = (other.numj,
                         other.nendo,other.nexo,
                         other.ncon,other.sigma,
                         other.jAA,other.jBB,
                         other.vardic,other.vdic,
                         other.modname,other.audic,
                         other.numjl,
                         other.nother)
                else:
                    intup = (other.numj,
                         other.nendo,other.nexo,
                         other.ncon,other.sigma,
                         other.jAA,other.jBB,
                         other.vardic,other.vdic,
                         other.modname,other.audic)
                other.modsolvers.forkleind = ForKleinD(intup)
################## 2ND-ORDER NON-LINEAR METHODS !!! ##################
#                if sum([nreg.search(x)!=None for x in other.txtpars.secs['vcvm'][0]]) == 0 and\
                if any([False if 'None' in x else True for x in secs['manualss'][0]]) and\
                   'numh' in dir(other):
                    # Open the MatKlein2D object
                    if 'nlsubsys' in dir(other):
                        intup = (other.numj,other.numh,
                             other.nendo,other.nexo,
                             other.ncon,other.sigma,
                             other.jAA,other.jBB,
                             other.vardic,other.vdic,
                             other.modname,other.audic,
                             other.numjl,other.numhl,
                             other.nother,sess1)
                    else:
                        intup = (other.numj,other.numh,
                             other.nendo,other.nexo,
                             other.ncon,other.sigma,
                             other.jAA,other.jBB,
                             other.vardic,other.vdic,
                             other.modname,other.audic,
                             sess1)
                    other.modsolvers.matklein2d = MatKlein2D(intup)
                    # Open the PyKlein2D object
                    intup = intup[:-1]
                    other.modsolvers.pyklein2d = PyKlein2D(intup)
    
            if 'jAA' in dir(other):
                other.mkeigv()

    def __getitem__(self,key):
        return self.wrapdic[key]
    def __repr__(self):
        return self.wrapdic.__repr__()
    def __str__(self):
        return self.wrapdic.__str__()
