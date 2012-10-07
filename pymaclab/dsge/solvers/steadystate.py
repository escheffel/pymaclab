'''
.. module:: steadystate
   :platform: Linux
   :synopsis: The steadystate module currently contains a shell class for sub-branching as well as two classes for calculating the
              steady state of DSGE models. One of them computes steady states using closed form information, while the other one makes
              use of a nonlinear root-finding algorithm supplied by the scipy optimize package.

.. moduleauthor:: Eric M. Scheffel <eric.scheffel@nottingham.edu.cn>


'''
import re
import numpy as np
import copy
from scipy import optimize
# Switch off runtime warnings here
import warnings
warnings.filterwarnings("ignore")

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
        rlog = re.compile(_rlog)
        rexp = re.compile(_rexp)
        manss={}
        locals().update(self.paramdic)
        globals().update(self.paramdic)
        for x in list_tmp1:
            str_tmp3 = x[:]
            while rlog.search(str_tmp3):
                str_tmp3 = re.sub(rlog,'np.log(',str_tmp3)
            while rexp.search(str_tmp3):
                str_tmp3 = re.sub(rexp,'np.exp(',str_tmp3)           
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
        # create a dic holding non _bar values in ssi
        tmp_dic = {}
        for keyo in self.ssi.keys():
            if '_bar' not in keyo: tmp_dic[keyo] = self.ssi[keyo]
        self.nonbardic = tmp_dic

    def solve(self,bounds_dic=None):

        # Turn the non-linear system into a representation
        # that is suitable for fsolve(invar) !

        ssi = self.ssi
        subdic = {}
        for y,z in zip(ssi.items(),range(len(ssi.items()))):
            subdic[y[0]] = (y[0],y[1],z)
            
        list_tmp1 = copy.deepcopy(self.ssm)
        for var in ssi.keys():
            _mreg = '(\+|\*|-|/|^|\()'+var
            mreg = re.compile(_mreg)
            _rlog = 'LOG\('
            _rexp = 'EXP\('
            rlog = re.compile(_rlog)
            rexp = re.compile(_rexp)
            for i1,line in enumerate(list_tmp1):
                # First strip away any whitespace
                while ' ' in list_tmp1[i1]:
                    list_tmp1[i1] = list_tmp1[i1].replace(' ','')
                while rlog.search(list_tmp1[i1]):
                    list_tmp1[i1] = re.sub(rlog,'np.log(',list_tmp1[i1])
                while rexp.search(list_tmp1[i1]):
                    list_tmp1[i1] = re.sub(rexp,'np.exp(',list_tmp1[i1])                 
                while mreg.search(list_tmp1[i1]):
                    ma = mreg.search(list_tmp1[i1])
                    matot = ma.group()
                    if len(ma.group(1))>0:
                        var = ma.group()[1:]
                        pos = ma.span()[0]+1
                    else:
                        var = ma.group()
                        pos = ma.span()[0]                      
                    poe = ma.span()[1]
                    list_tmp1[i1] = list_tmp1[i1][:pos]+'invar['+str(subdic[var][2])+']'+list_tmp1[i1][poe:]
        func_repr = list_tmp1
        self.subdic = copy.deepcopy(subdic)
        self.ssm_alt = copy.deepcopy(func_repr)

        # Define the function to be handed over
        # to fsolve
        if bounds_dic == None:
            def func(invar):
                locals().update(self.paramdic)
                locals().update(self.nonbardic)
                fdot = np.zeros((len(func_repr)),float)
                for i1,x in enumerate(func_repr):
                    fdot[i1] = eval(x)
                return fdot
        else:
            def func(invar):
                locals().update(self.paramdic)
                locals().update(self.nonbardic)
                fdot = np.zeros((len(func_repr)),float)
                for i1,x in enumerate(func_repr):
                    fdot[i1] = eval(x)
                return np.sum(np.abs(fdot))
            

        # Define the initial values and
        # start the non-linear solver
        inlist = []
        for x in subdic.values():
            inlist.append([x[2],x[1],x[0]])
        inlist.sort()
        init_val = [float(x[1]) for x in inlist]

        # Prepare the case of constrained optimisation if needed
        blist = []
        if bounds_dic != None:
            for elem in inlist:
                if elem[2] in bounds_dic.keys(): blist.append(bounds_dic[elem[2]])
                else: blist.append((None,None))

        # Accomodate unconstrained and constrained optimisation
        if bounds_dic == None:
            outobj = optimize.root(func,init_val,method='hybr')
            output = outobj.x
            mesg = outobj.message
            ier = outobj.status
        else:
            outobj = optimize.minimize(func,init_val,method='L-BFGS-B',bounds=blist)
            output = outobj.x
            mesg = outobj.message
            ier = outobj.status
        # Attach the outputs of the solver as attributes
        self.fsout={}
        for x,y in zip(output,inlist):
            self.fsout[y[2]] = x
        self.ier = ier
        self.mesg = mesg
        self.output = output
