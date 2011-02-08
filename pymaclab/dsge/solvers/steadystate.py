"""
THE STEADY STATE SOLVER CLASS AND ITS SUBCLASSES
"""
import re
import numpy as np
import copy
from scipy import optimize

class SSsolvers(object):
    def __init__(self):
        pass

class ManualSteadyState(SSsolvers):
    """
    Solve for the steady state manually
    """
    def __init__(self, intup):
        self.manss_sys = intup[0]
        self.paramdic = intup[1]

    def solve(self):
        list_tmp1 = copy.deepcopy(self.manss_sys)
        # Create manual (closed-form) steady state dictionary
        _rlog = 'LOG\('
        _rexp = 'EXP\('
        _fsolve = 'ROOT\((?P<nlexp>.*),\s*(?P<vari>.*)\s*=\s*(?P<init>.*)\s*,\s*fail\s*=\s*(?P<fval>.*)\);'
        rlog = re.compile(_rlog)
        rexp = re.compile(_rexp)
        fexp = re.compile(_fsolve)
        manss={}
        locals().update(self.paramdic)
        globals().update(self.paramdic)
        for x in list_tmp1:
            str_tmp3 = x[:]
            while rlog.search(str_tmp3):
                str_tmp3 = re.sub(rlog,'np.log(',str_tmp3)
            while rexp.search(str_tmp3):
                str_tmp3 = re.sub(rexp,'np.exp(',str_tmp3)           

            # Do fsolve root finding, if ROOT detected
            #NOTE: use ROOT() to solve non-linear FOCs (see modfiles/max2.txt)
            if fexp.search(str_tmp3):
                asvari = str_tmp3.split('=')[0].strip()
                ma = fexp.search(str_tmp3)
                nlexp = ma.group('nlexp')
                vari = ma.group('vari')
                init = np.float(ma.group('init'))
                fval = np.float(ma.group('init'))
                solu,infodict,ier,mesg = \
                    optimize.fsolve(eval('lambda '+vari+':'+nlexp),
                                        init,full_output=1)
                if ier == 1:
                    manss[asvari] = solu
                    locals()[asvari] = solu
                    globals().update(manss)
                else:
                    manss[asvari] = fval
                    locals()[asvari] = fval
                    globals().update(manss)
                continue
            str_tmp = str_tmp3.split(';')[0]
            list_tmp = str_tmp.split('=')
            str_tmp1 = list_tmp[0].strip()
            str_tmp2 = list_tmp[1].strip()
            manss[str_tmp1] = eval(str_tmp2)
            locals()[str_tmp1] = eval(str_tmp2)
            globals().update(manss)
        self.sstate = manss

class Fsolve(SSsolvers):
    def __init__(self,intup):
        self.ssm = intup[0]
        self.ssi = intup[1]
        self.paramdic = intup[2]

    def solve(self):

        # Turn the non-linear system into a representation
        # that is suitable for fsolve(invar) !

        ssi = self.ssi
        subdic = {}
        for y,z in zip(ssi.items(),range(len(ssi.items()))):
            subdic[y[0]] = (y[0],y[1],z)
        list_tmp1 = copy.deepcopy(self.ssm)
        for var in self.ssi.keys():
            _mreg = '(\+|\*|-|/|^|[(])'+var
            mreg = re.compile(_mreg)
            for i1,line in enumerate(list_tmp1):
                while mreg.search(list_tmp1[i1]):
                    ma = mreg.search(list_tmp1[i1])
                    if len(ma.group(1))>0:
                        var = ma.group()[1:]
                        pos = ma.span()[0]+1
                    else:
                        var = ma.group()
                        pos = ma.span()[0]                      
                    poe = ma.span()[1]
                    list_tmp1[i1] = list_tmp1[i1][:pos]+'invar['+str(subdic[var][2])+']'+list_tmp1[i1][poe:]

        func_repr = list_tmp1

        # Define the function to be handed over
        # to fsolve
        def func(invar):
            locals().update(self.paramdic)
            fdot = np.zeros((len(func_repr)),float)
            i1=0
            for x in func_repr:
                fdot[i1] = eval(x)
                i1 = i1 + 1
            return fdot

        # Define the initial values and
        # start the non-linear solver
        inlist = []
        for x in subdic.values():
            inlist.append([x[2],x[1],x[0]])
        inlist.sort()
        init_val = [x[1] for x in inlist]
        (output,infodict,ier,mesg) = optimize.fsolve(func,init_val,
                                                full_output=1)
        # Attach the outputs of the solver as attributes
        self.fsout={}
        for x,y in zip(output,inlist):
            self.fsout[y[2]] = x
        self.infodict = infodict
        self.ier = ier
        self.mesg = mesg
