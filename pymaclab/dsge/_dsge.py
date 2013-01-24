'''
.. module:: _dsge
   :platform: Linux
   :synopsis: The _dsge module contains the the most important class for doing work with DSGE models.
              In particular it supplies the DSGEmodel class containing most of the functionality of
              DSGE model instances.

.. moduleauthor:: Eric M. Scheffel <eric.scheffel@nottingham.edu.cn>


'''
##############################################################################
##                       IMPORTING OF EXTERNAL LIBRARIES                    ##
##############################################################################
#from __future__ import division
import copy
from copy import deepcopy
import os
import subprocess
import sys
import re
from pymaclab.dsge.helpers._helpers import now_is
import numpy as np
from numpy import matlib as mat
from scipy.linalg import eig as scipyeig
from scipy import optimize
from tempfile import mkstemp
import time
import datetime
import glob
# Import Sympycore, now comes supplied and installed with pymaclab
import sympycore
# Import PP, now comes supplied and installed with pymaclab
import pp
##############################################################################
##                  IMPORTING OF EXTERNAL LIBRARIES DONE!                   ##
##############################################################################


##############################################################################
##                    IMPORTING REFACTORED CLASSES/METHODS                  ##
##############################################################################
# Import of refactored filtering routines
from pymaclab.filters._hpfilter import hpfilter
from pymaclab.filters._bkfilter import bkfilter

# Import of refactored statistics Classes
from pymaclab.stats.var import VAR

# Import refactored DSGE model initialisor class and the dynarepp flag
from pymaclab.dsge.inits._inits import Inits
import pymaclab.dsge.inits._set_flags
from pymaclab.dsge.inits._set_flags import dynare_flag

# Import refactored additional helpful objects
from pymaclab.dsge.helpers._paths import *
from pymaclab.dsge.helpers._greek_alphab import greek_alph

# Make sure the jobserver has done his global jobs each time you instantiate a new instance
global jobserver,ppservers
from pymaclab.dsge.inits._inits import jobserver,ppservers
##############################################################################
##               IMPORTING REFACTORED CLASSES/METHODS DONE!                 ##
##############################################################################
    


##################THE DSGE MODEL CLASS (WORKS)#####################
class DSGEmodel(object):
    '''
    This is the macrolab DSGEmodel class. It is the main class of the packages and instantiates
    DSGE model instances which possess many features such as model file parsers, solvers, etc. The __init__
    function first called mostly attaches the passed arguments to the DSGE model instance in form of private _X
    data fields.
    
    .. note::
    
       Notice that the various init method, i.e. init1, init1a, etc. are not called here, but they are called externally in the
       pymaclab package when instantiating a new DSGE model using pymaclab.newMOD().
    
    :param ffile:             The absolute path to the PyMacLab DSGE model file to be parsed on instantiation
    :type ffile:              str
    :param dbase:             A database with time series data for estimation purposes. Not implemented at the moment
    :type standard:           unknown
    :param initlev:           Takes values 0,1,2. Determines how deep the instantiation cascades through all methods.
                              If 0 then the file is parsed and the instance is only *prepared* for steady state solving.
                              If 1 then the file is parsed, the SS is automatically computed and instance is *prepared*
                              for dynamic solving. If 2 then the instance is solved all the way.
    :type initlev:            int
    :param mesg:              If True then lots of diagnostics are printed to the screen during instantiation.
    :type mesg:               bool
    :param ncpus:             The number of CPU core to be employed, defaults to 1. But 'auto' can also be used for detection
    :type ncpus:              int|str
    :param mk_hessian:        Whether the Hessian should be computed as this is expensive.
    :type mk_hessian:         bool
    :param use_focs:          Should the FOCs be used directly to look for the steady state? Must use list|tuple to pick equations.
    :type use_focs:           tuple
    :param ssidic:            A Python ssidic with the initial starting values for solving for the SS numerically
    :type ssidic:             dic
    :param sstate:            A Python ssidic with the externally computed steady state values of the model
    :type sstate:             dic
        
    :return self:             *(dsge_inst)* - A populated DSGE model instance with fields and methods

    '''
    class __Derivatives(object):
        '''
        Empty inner shell class for saving derivatives class
        '''
        pass

    # Initializes the absolute basics, errors difficult to occur
    def __init__(self,ffile=None,dbase=None,initlev=2,mesg=False,ncpus='auto',mk_hessian=True,\
                 use_focs=False,ssidic=None,sstate=None,vtiming={'exo':[-1,0],'endo':[-1,0],'con':[0,1]}):
        jobserver.wait()
        if sstate != None:
            self._sstate = deepcopy(sstate)
        self._initlev = initlev #TODO: remove this as an option
        self._vtiming = vtiming
        self._mesg = mesg
        self._ncpus = ncpus
        self._mk_hessian = mk_hessian
        self._use_focs = use_focs
        self._ssidic = ssidic
        if self._use_focs and self._ssidic != None:
            self.ssidic = deepcopy(self._ssidic)
        # Set no author
        self.setauthor()
        # Open the switches dic and initiate
        #NOTE this is _not_ pythonic
        self.switches = {}
        self.switches['errocc'] = '0'
        self.switches['ss_suc'] = ['0','0']
        # Attach some model attributes
        self.modfile = ffile
        self.dbase = dbase
        # Open the derivates branch using inner class to save derivatives there
        self.derivatives = DSGEmodel.__Derivatives()         
        # Open the inits branch at the end
        self.inits = Inits(other=self)

            
    def find_rss(self,mesg=False,rootm='hybr',scale=0.0):
        '''
        The is a method which can be called to find the risky steady state
           
        :param self:     The DSGE model instance itself.
        :type self:      dsge_inst

        '''
         
        varbar = []
        nexo = self.nexo
        for elem in self.vardic['endo']['var']:
            varbar.append(elem[0].split('(')[0]+'_bar')
        for elem in self.vardic['con']['var']:
            varbar.append(elem[0].split('(')[0]+'_bar')
        tmp_dic = {}
        tmp_dic.update(self.paramdic)
        tmp_dic.update(self.sstate)
        sstate_li = []
        sstate = {}
        for elem in varbar:
            if elem in tmp_dic.keys():
                sstate[elem] = tmp_dic[elem]
                sstate_li.append(tmp_dic[elem])
        rsstate_li = deepcopy(sstate_li)
        rsstate = deepcopy(sstate)

        clone = copy.deepcopy(self)
        clone._mesg = False
        
        # Define the function to be handed over
        # to fsolve
        def func(invar):
            '''
            ####### Put in a trap in case number turn negative ########
            negli = [True if x < 0.0 else False for x in invar]
            if any(negli):
                for i1,elem in enumerate(negli):
                    invar[i1] = self.sstate[varbar[i1]]+0.01*self.sstate[varbar[i1]]
            ###########################################################
            '''
            # Update ss with passed values
            for i1,elem in enumerate(varbar):            
                clone.sstate[elem] = invar[i1]
            # Update the model's derivatives, but only numerically
            clone.init5(update=True)
            # Solve the 2nd-order accurate solution
            clone.modsolvers.pyklein2d.solve()
            # Get retval into right shape
            KX = clone.modsolvers.pyklein2d.KX
            retval = [float(x) for x in KX[nexo:]]
            KY = clone.modsolvers.pyklein2d.KY
            retval+=[float(x) for x in KY]
            retval = np.array(retval)
            if mesg:
                print invar
                print '------------------------------------'
            return retval

        
        # Define the initial values and
        # start the non-linear solver
        init_val = np.array([float(self.sstate[varbar[i1]])+\
                             (scale/100.0)*float(self.sstate[varbar[i1]]) for i1 in range(len(varbar))])
        outobj = optimize.root(func,init_val,method=rootm)
        rss_funcval1 = outobj.fun
        output = outobj.x
        mesg = outobj.message
        ier = outobj.status
        
        self.rss_funcval1 = rss_funcval1
                
        # Check the difference and do bounded minimisation (root-finding)
        diff_dic = {}
        bounds_dic = {}
        for i1,keyo in enumerate(varbar):
            if keyo in self.sstate.keys():
                if output[i1] != self.sstate[keyo]:
                    diff_dic[keyo] = output[i1]
                    bounds_dic[keyo] = (output[i1],output[i1])
        if bounds_dic != {}:
            # Also add the non-bar SS values to bounds_dic, before doing constrained minimisation
            for keyo in self.sstate.keys():
                if "_bar" not in keyo: bounds_dic[keyo] = (self.sstate[keyo],self.sstate[keyo])
            # Do constrained root finding here
            if 'fsolve' in dir(clone.sssolvers):
                clone.sssolvers.fsolve.solve(bounds_dic=bounds_dic)
            # Also need to make sure residually computed SS get updated accordingly
            if 'manss' in dir(clone.sssolvers) and 'fsolve' in dir(clone.sssolvers):
                    clone.sssolvers.manss.paramdic.update(clone.sssolvers.fsolve.fsout)
                    clone.sssolvers.manss.solve()

            # Run loop one last time to get rss_funval2
            if 'fsolve' in dir(clone.sssolvers):
                clone.sstate.update(clone.sssolvers.fsolve.fsout)
            if 'manss' in dir(clone.sssolvers):
                clone.sstate.update(clone.sssolvers.manss.sstate)
            clone.init5(update=True)
            clone.modsolvers.pyklein2d.solve()
            retval = []
            # Get retval into right shape
            KX = clone.modsolvers.pyklein2d.KX
            retval = [float(x) for x in KX]
            KY = clone.modsolvers.pyklein2d.KY
            for elem in KY:
                retval.append(float(elem))
            retval = np.array(retval)
            self.rss_funcval2 = retval
             
        # Now attach final results to instance
        self.rsstate = deepcopy(self.sstate)
        self.rparamdic = deepcopy(self.paramdic)
        for i1,elem in enumerate(varbar):
            if elem in self.rsstate.keys(): self.rsstate[elem] = output[i1]
            if elem in self.rparamdic.keys(): self.rparamdic[elem] = output[i1]
        if bounds_dic != {}:
            if 'fsolve' in dir(clone.sssolvers):
                self.rsstate.update(clone.sssolvers.fsolve.fsout)
            if 'fsolve' in dir(clone.sssolvers) and 'manss' in dir(clone.sssolvers):
                self.rsstate.update(clone.sssolvers.manss.sstate)

            
        

    # html info of model opened with webbrowser
    def info(self):
        '''
        A convience method for collecting and displaying the model's properties in a browser
        using html language.
        '''
        tmplist = glob.glob('tempz*.html')
        for x in tmplist:
            os.remove(x)
        modname = self.mod_name
        secs = self.txtpars.secs
        direc = os.getcwd()
        fd,fpath = mkstemp(prefix='tempz',suffix='.html',dir=direc)
        htmlfile = os.fdopen(fd,'w+b')
        htmlfile.write('<HTML><BODY BGCOLOR="white">\n')
        htmlfile.write('<H2>%s</H2>\n'%modname)
        htmlfile.write('\n\n')
        for x1 in secs['info'][1]:
            htmlfile.write('<P>'+x1+'\n')
        htmlfile.write('<P>'+'<H4>Model Parameters</H4>\n')
        for x1 in secs['params'][1]:
            htmlfile.write('<P>'+x1+'\n')
        htmlfile.write('<P>'+'<H4>Nonlinear First-Order Conditions</H4>\n')
        for x1 in secs['focs'][1]:
            if x1[0] == '[':
                htmlfile.write('<P>'+x1+'\n')
            elif x1[0] != '[':
                htmlfile.write('<P>'+4*'&nbsp;'+ x1+'\n')
        htmlfile.write('<P>'+'<H4>Log-Linearized Model Equations</H4>\n')
        for x1 in secs['modeq'][1]:
            if x1[0] == '[':
                htmlfile.write('<P>'+x1+'\n')
            elif x1[0] != '[':
                htmlfile.write('<P>'+4*'&nbsp;'+ x1+'\n')
        htmlfile.write('</BODY></HTML>\n')
        htmlfile.write('<P>'+'<H2>Solution Section</H2>\n')
        htmlfile.write('<P>'+'<H4>Uhlig Toolkit Output</H4>\n')
        htmlfile.close()
        cmd = browserpath+' '+fpath
        os.system(cmd)
        tmplist = glob.glob('tempz*.html')
        for x in tmplist:
            os.remove(x)
        return 'Model Website opened!'
    
    def setauthor(self,author=None):
        '''
        Convience method for setting the model's author's name.
        '''
        if author == None:
            self.author = 'No author'
        else:
            self.author = author

    def pdf(self):
        '''
        This will pdflatex the model's tex file and then view it in the system's pdf viewer.
        '''
        modfile = self.modfile
        modfname = modfile.split('.')[0]
        if modfname+'.pdf' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.pdf'))
        if modfname+'.log' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.log'))
        if modfname+'.aux' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.aux'))
        if modfname+'.dvi' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.dvi'))
        if modfname+'.log' in os.listdir(os.getcwd()):
            os.remove(os.path.join(os.getcwd(),modfname+'.log'))

        # Does the model tex file even exist?
        if modfname+'.tex' not in os.listdir(modfpath):
            print 'Error: The model tex file does not exist!'
            print 'Use model.texed() to create a new one!'
            return

        # First check for Chktex syntax errors
        if 'texer.log' in os.listdir(os.getcwd()):
            os.remove(os.path.join(os.getcwd(),'texer.log'))
        args = '--quiet -v2 -n3 -n25 -n12 -n35'
        cmd = 'chktex '+args+' '+os.path.join(modfpath,modfname+'.tex'+' > texer.log')
        os.system(cmd)
        file = open(os.path.join(os.getcwd(),'texer.log'),'rU')
        errlist = file.readlines()
        if len(errlist)== 0:
            file.close()
            os.remove(os.path.join(os.getcwd(),'texer.log'))
        else:
            print 'There were errors in the model tex file! Abort!\n'

            for x in errlist:
                print x
            file.close()
            os.remove(os.path.join(os.getcwd(),'texer.log'))
            return

        # PDFLatex the tex file and open
        cmd = 'pdflatex '+'-output-directory='+modfpath+' '+os.path.join(modfpath,modfname+'.tex'+' > out.log')
        os.system(cmd)
        cmd2 = pdfpath + ' ' + os.path.join(modfpath,modfname+'.pdf')
        os.system(cmd2)
        if modfname+'.pdf' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.pdf'))
        if modfname+'.log' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.log'))
        if modfname+'.aux' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.aux'))
        if modfname+'.dvi' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.dvi'))
        if modfname+'.log' in os.listdir(os.getcwd()):
            os.remove(os.path.join(os.getcwd(),modfname+'.log'))

    def txted(self):
        '''
        A convience method for launching an editor editing the model's associated
        model txt file. If the file gets saved (even if unaltered) the model is
        re-initialized.
        '''
        modfile = self.modfile
        modfname = modfile.split('.')[0]
        if modfname+'.log' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.log'))
        if modfname+'.log' in os.listdir(os.getcwd()):
            os.remove(os.path.join(os.getcwd(),modfname+'.log'))
        cmd = txtedpath+' '+os.path.join(modfpath,self.modfile+' > out.log')
        timestamp0 = os.stat(os.path.join(modfpath,self.modfile))[8]
        os.system(cmd)
        timestamp1 = os.stat(os.path.join(modfpath,self.modfile))[8]
        if timestamp0 != timestamp1:
            self.__init__(self.modfile,self.dbase)
            self.init2()
        if modfname+'.log' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.log'))
        if modfname+'.log' in os.listdir(os.getcwd()):
            os.remove(os.path.join(os.getcwd(),modfname+'.log'))

    def texed(self):
        '''
        A convenience method which allows users to launch the model's tex file in an
        editor specified in the configuration settings of pymaclab.
        '''
        modfile = self.modfile
        modfname = modfile.split('.')[0]
        if modfname+'.log' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.log'))
        if modfname+'.log' in os.listdir(os.getcwd()):
            os.remove(os.path.join(os.getcwd(),modfname+'.log'))
        if modfname+'.tex' not in os.listdir(modfpath):
            file = open(os.path.join(modfpath,modfname+'.tex'),'wb')
            file.write('\\documentclass[a4paper,11pt]{article}\n')
            file.write('\\title{%s}\n'% self.mod_name)
            file.write('\\author{%s}\n'% self.author)
            file.write('\\begin{document}\n')
            file.write('\\maketitle\n\n')
            file.write('Start writing your model in Latex here!')
            file.write(8*'\n')
            file.write('\\end{document}\n')
            file.close()
        cmd = texedpath+' '+os.path.join(modfpath,self.modfile[:-3]+'tex'+' > out.log')
        os.system(cmd)
        if modfname+'.log' in os.listdir(modfpath):
            os.remove(os.path.join(modfpath,modfname+'.log'))
        if modfname+'.log' in os.listdir(os.getcwd()):
            os.remove(os.path.join(os.getcwd(),modfname+'.log'))

    def deltex(self):
        '''
        A simple method which allows you to delete the current tex file associated with this model.
        Will ask for confirmation though, so cannot be used in a fire-and-forget batch file.
        '''
        answ = 0
        while answ not in ['Y','y','N','n']:
            answ = raw_input("REALLY delete this model's tex file ? (Y/N)")
        if answ in ['Y','y']:
            modfile = self.modfile
            modfname = modfile.split('.')[0]
            if modfname+'.tex' in os.listdir(modfpath):
                os.remove(os.path.join(modfpath,modfname+'.tex'))
        elif answ in ['N','n']:
            return

    def ccv(self,instring=None):
        '''
        This is just a convenience method which allows the current preferred dynamic
        solution method to be stored (attached) here and made accessible without having
        to traverse too deeply into the model instance tree structure.
        '''
        if instring == None:
            return 'Error: Your have not specified a current view name!'
        elif instring not in dir(self.modsolvers):
            return 'Error: Your chosen view does not exist for this model!'
        else:
            self.cv = eval('self.modsolvers.'+instring)

    def vreg(self,paratuple=(None,'all','0'),cinstring='',iter=False,info='min'):
        '''
        The regex function for variable detection out of strings
        '''
        vars = {}
        vars['con'] = []
        vars['endo'] = []
        vars['exo'] = []
        vars['other'] = []
        vars['iid'] = []
        vars['state'] = []
        vars['all'] = []
        instring = deepcopy(cinstring)

        for x in self.vardic['con']['var']:
            vars['con'].append(x[0].split('(')[0])
        for x in self.vardic['endo']['var']:
            vars['endo'].append(x[0].split('(')[0])
        for x in self.vardic['exo']['var']:
            vars['exo'].append(x[0].split('(')[0])
            vars['iid'].append(x[2].split('(')[0])
        for x in self.vardic['other']['var']:
            if '|' in x[0]:
                vars['other'].append(x[0].split('|')[1].split('(')[0])
            else:
                vars['other'].append(x[0].split('(')[0])
        vars['state'] = vars['endo'] + vars['exo']
        vars['all'] = vars['state'] + vars['con'] + vars['iid'] + vars['other']

        # Create regex part for expectation info set
        if paratuple[0] != None and '{' not in paratuple[0]:
            if paratuple[0] == '0':
                iti = '(E\(t(?P<itime>.{0,0})\)\|)'
            elif paratuple[0][0] == '-':
                iti = '(E\(t(?P<itime>-'+paratuple[0][1:]+')\)\|)'
            else:
                iti = '(E\(t(?P<itime>\+'+paratuple[0]+')\)\|)'
        elif paratuple[0] != None and '{' in paratuple[0]:
            if '|' not in paratuple[0]:
                indx1 = paratuple[0].split(',')[0][1:]
                indx2 = paratuple[0].split(',')[1][:-1]
                indx2 = str(eval(indx2+'+1'))
            elif '|' in paratuple[0]:
                strtmp = paratuple[0].split('|')[0]
                indx1 = strtmp.split(',')[0][1:]
                indx2 = strtmp.split(',')[1][:-1]
                indx2 = str(eval(indx2+'+1'))
            itrange = eval('range('+indx1+','+indx2+',1)')
            zswitch = 0
            if 0 in itrange:
                zswitch = 1
                itrange.pop(itrange.index(0))       
            str_tmp=''
            for x in itrange:
                if str(x)[0] == '-':
                    str_tmp = str_tmp + '|'+str(x)
                elif str(x) == '0':
                    str_tmp = str_tmp + '|'+str(x)
                else:
                    str_tmp = str_tmp + '|\+'+str(x)
            str_tmp = str_tmp[1:]+zswitch*'|.{0,0}'
            orswitch = False
            if '|' in paratuple[0]:
                orswitch = True
            iti = '(E\(t(?P<itime>'+str_tmp[:]+')\)\|)'+(1-orswitch)*'{1,1}'+orswitch*'{0,1}'
        elif paratuple[0] == None:
            iti = ''


        # Create regex part for variable type
        if '|' not in paratuple[1]:
            varp = ''
            for x in vars[paratuple[1]]:
                varp = varp + '|'+x
        elif '|' in paratuple[1]:
            typeli = paratuple[1].split('|')
            varp = ''
            for x in typeli:
                for y in vars[x]:
                    varp = varp + '|'+y
        varp = varp[1:]
        varp = '(?P<var>'+varp+')'


        # Create the time part
        str_tmp = ''
        if '{' not in paratuple[2]:
            if paratuple[2][0] == '-':
                ti = '\(t(?P<time>'+paratuple[2]+')\)'
            elif paratuple[2] == '0':
                ti = '\(t(?P<time>.{0,0})\)'
            else:
                ti = '\(t(?P<time>\+'+paratuple[2]+')\)'
        elif '{' in paratuple[2]:
            indx1 = paratuple[2].split(',')[0][1:]
            indx2 = paratuple[2].split(',')[1][:-1]
            indx2 = str(eval(indx2+'+1'))
            trange = eval('range('+indx1+','+indx2+',1)')
            zswitch = False
            if 0 in trange:
                zswitch = True
                trange.pop(trange.index(0))
            str_tmp=''
            for x in trange:
                if str(x)[0] == '-':
                    str_tmp = str_tmp + '|'+str(x)
                elif str(x) == '0':
                    str_tmp = str_tmp + '|'+str(x)
                else:
                    str_tmp = str_tmp + '|\+'+str(x)
            str_tmp = str_tmp[1:]+zswitch*'|.{0,0}'
            ti = '(\(t(?P<time>'+str_tmp[:]+')\))'

        reva = re.compile('(?P<pre>^[^a-zA-Z^_]{0,1}|[^a-zA-Z^_]{1,1})('+iti+varp+ti+')')

        if reva.search(instring) and iter == False:
            ma = reva.search(instring)
            if ma.group('pre') != None:
                pre = ma.group('pre')
                prel = len(pre)
            else:
                prel = 0
            alcap = ma.group()[prel:]
            var = ma.group('var')
            time = ma.group('time')
            if time != '':
                time = ma.group('time')
            else:
                time = '0'
            span = (ma.span()[0]+prel,ma.span()[1])
            vtype = (var in vars['endo'])*'endo'+\
                  (var in vars['exo'])*'exo'+\
                  (var in vars['con'])*'con'+\
                  (var in vars['iid'])*'iid'+\
                  (var in vars['other'])*'other'
            if vtype == 'endo':
                posx = [x[0].split('(')[0] for x in self.vardic['endo']['var']].index(var)
            elif vtype == 'exo':
                posx = [x[0].split('(')[0] for x in self.vardic['exo']['var']].index(var)
            elif vtype == 'con':
                posx = [x[0].split('(')[0] for x in self.vardic['con']['var']].index(var)
            elif vtype == 'iid':
                posx = [x[2].split('(')[0] for x in self.vardic['exo']['var']].index(var)
            elif vtype == 'other':
                posx = [x[0] for x in self.vardic['other']['var']].index(alcap)

            if ma.groupdict().has_key('itime'):
                itime = ma.group('itime')
                if itime != '':
                    itime = ma.group('itime')
                else:
                    itime = '0'
            else:
                itime = None

            if info == 'max':
                return (alcap,(vtype,posx),(itime,var,time),span)
            elif info == 'min':
                return (alcap,(vtype,posx),(itime,var,time))
        elif reva.search(instring) and iter == True:
            itobj = reva.finditer(instring)
            relist=[]
            for ma in itobj:
                if ma.group('pre') != None:
                    pre = ma.group('pre')
                    prel = len(pre)
                else:
                    prel = 0
                alcap = ma.group()[prel:]
                var = ma.group('var')
                time = ma.group('time')
                if time != '':
                    time = ma.group('time')
                else:
                    time = '0'
                span = (ma.span()[0]+prel,ma.span()[1])
                vtype = (var in vars['endo'])*'endo'+\
                      (var in vars['exo'])*'exo'+\
                      (var in vars['con'])*'con'+\
                      (var in vars['iid'])*'iid'+\
                      (var in vars['other'])*'other'
                if vtype == 'endo':
                    posx = [x[0].split('(')[0] for x in self.vardic['endo']['var']].index(var)
                elif vtype == 'exo':
                    posx = [x[0].split('(')[0] for x in self.vardic['exo']['var']].index(var)
                elif vtype == 'con':
                    posx = [x[0].split('(')[0] for x in self.vardic['con']['var']].index(var)
                elif vtype == 'iid':
                    posx = [x[2].split('(')[0] for x in self.vardic['exo']['var']].index(var)
                elif vtype == 'other':
                    posx = [x[0] for x in self.vardic['other']['var']].index(alcap)
                if ma.groupdict().has_key('itime'):
                    itime = ma.group('itime')
                    if itime != '':
                        itime = ma.group('itime')
                    else:
                        itime = '0'
                else:
                    itime = None

                relist.append((alcap,(vtype,posx),(itime,var,time),span))
            if info == 'min':
                relist2 = [(x[0],x[1],x[2]) for x in relist]
                for x in relist2:
                    while relist2.count(x) > 1:
                        indx = relist2.index(x)
                        relist2.pop(indx)
                return relist2
            elif info == 'max':
                return relist
        else:
            return False
###########ANALYTIC AND NUMERICAL JACOBIAN AND HESSIAN METHODS############
    def def_differ_periods(self):
        '''
        A convenience method which uses the vtiming attribute in order to create a number of lists of
        variables with their time definitions (subscripts) for which derivates have to be taken using the
        non-linear system of FOCs.
        
        :param self: object instance
        :type self: dsge_inst
        
        :return self.jAA:  *(arr2d)* - attaches numerical AA matrix used in Forkleind solution method
        :return self.jBB:  *(arr2d)* - attaches numerical BB matrix used in Forkleind solution method
        '''
        vtiming = deepcopy(self.vtiming)
        # Timing assumptions, first for exogenous variables
        # For past
        if vtiming['exo'][0] < 0:
            exo_0 = [x[0].split('(')[0]+'(t-'+str(abs(vtiming['exo'][0]))+')' for x in self.vardic['exo']['var']]
        elif vtiming['exo'][0] == 0:
            exo_0 = [x[0].split('(')[0]+'(t)' for x in self.vardic['exo']['var']]
        elif vtiming['exo'][0] > 0:
            exo_0 = ['E(t)|'+x[0].split('(')[0]+'(t+'+str(vtiming['exo'][0])+')' for x in self.vardic['exo']['var']]
        # For future
        if vtiming['exo'][1] < 0:
            exo_1 = [x[0].split('(')[0]+'(t-'+str(abs(vtiming['exo'][1]))+')' for x in self.vardic['exo']['var']]
        elif vtiming['exo'][1] == 0:
            exo_1 = [x[0].split('(')[0]+'(t)' for x in self.vardic['exo']['var']]
        elif vtiming['exo'][1] > 0:
            exo_1 = ['E(t)|'+x[0].split('(')[0]+'(t+'+str(vtiming['exo'][1])+')' for x in self.vardic['exo']['var']]
            
        # Timing assumptions, endogenous variables
        # For past
        if vtiming['endo'][0] < 0:
            endo_0 = [x[0].split('(')[0]+'(t-'+str(abs(vtiming['endo'][0]))+')' for x in self.vardic['endo']['var']]
        elif vtiming['endo'][0] == 0:
            endo_0 = [x[0].split('(')[0]+'(t)' for x in self.vardic['endo']['var']]
        elif vtiming['endo'][0] > 0:
            endo_0 = [x[0].split('(')[0]+'(t+'+str(vtiming['endo'][0])+')' for x in self.vardic['endo']['var']]
        # For future, BE CAREFUL, no expectations term on variables with (t+1) for endo
        if vtiming['endo'][1] < 0:
            endo_1 = [x[0].split('(')[0]+'(t-'+str(abs(vtiming['endo'][1]))+')' for x in self.vardic['endo']['var']]
        elif vtiming['endo'][1] == 0:
            endo_1 = [x[0].split('(')[0]+'(t)' for x in self.vardic['endo']['var']]
        elif vtiming['endo'][1] > 0:
            endo_1 = [x[0].split('(')[0]+'(t+'+str(vtiming['endo'][1])+')' for x in self.vardic['endo']['var']]
            
        # Timing assumptions, control variables
        # For past
        if vtiming['con'][0] < 0:
            con_0 = [x[0].split('(')[0]+'(t-'+str(abs(vtiming['con'][0]))+')' for x in self.vardic['con']['var']]
        elif vtiming['con'][0] == 0:
            con_0 = [x[0].split('(')[0]+'(t)' for x in self.vardic['con']['var']]
        if vtiming['con'][0] > 0:
            con_0 = ['E(t)|'+x[0].split('(')[0]+'(t+'+str(vtiming['con'][0])+')' for x in self.vardic['con']['var']]        
        # For future
        if vtiming['con'][1] < 0:
            con_1 = [x[0].split('(')[0]+'(t-'+str(abs(vtiming['con'][1]))+')' for x in self.vardic['con']['var']]
        elif vtiming['con'][1] == 0:
            con_1 = [x[0].split('(')[0]+'(t)' for x in self.vardic['con']['var']]
        elif vtiming['con'][1] > 0:
            con_1 = ['E(t)|'+x[0].split('(')[0]+'(t+'+str(vtiming['con'][1])+')' for x in self.vardic['con']['var']]
        # Also attach the final differvardic to the model for inspection
        differvardic = {}
        differvardic['exo_0'] = exo_0
        differvardic['exo_1'] = exo_1
        differvardic['endo_0'] = endo_0
        differvardic['endo_1'] = endo_1
        differvardic['con_0'] = con_0
        differvardic['con_1'] = con_1
        self.differvardic = deepcopy(differvardic)
        return exo_0,exo_1,endo_0,endo_1,con_0,con_1

    def mkjahe(self):
        '''
        An unparallelized method using native Python and Sympycore in oder
        to calculate the numerical and analytical Jacobian and Hessian of the model.
        
        :param self: object instance
        :type self: dsge_inst
        
        :return self.numj: *(arr2d)* - attaches numerical model Jacobian to instance
        :return self.jdic: *(dic)* - attaches the analytical model Jacobian to instance
        :return self.numh: *(arr2d)* - attaches numerical model Hessian to instance
        :return self.hdic: *(dic)* - attaches the analytical 3D Hessian to instance
        :return self.jAA:  *(arr2d)* - attaches numerical AA matrix used in Forkleind solution method
        :return self.jBB:  *(arr2d)* - attaches numerical BB matrix used in Forkleind solution method
        '''
        mk_hessian = self._mk_hessian
               
        exo_0,exo_1,endo_0,endo_1,con_0,con_1 = self.def_differ_periods()

        inlist = exo_1+endo_1+con_1+exo_0+endo_0+con_0
        intup=tuple(inlist)

        patup = ('{-10,10}|None','endo|con|exo|other','{-10,10}')
        jcols = len(intup)
        jrows = len(self.nlsys_list)
        nlsys = deepcopy(self.nlsys_list)

        if 'nlsubsys' in dir(self):
            lsubs = len(self.nlsubsys)
            jrows = jrows + lsubs
            nlsys = nlsys + deepcopy(self.nlsubsys)
            
        # Create and alternative list of variables for substitution purposes
        subsoli = []
        for lino in nlsys:
            list1 = self.vreg(patup,lino,True,'max')
            list1 = [x[0] for x in list1]
            for varo in list1:
                if varo not in subsoli: subsoli.append(varo)
        intup = deepcopy(subsoli)

        # Create substitution var list and dictionary
        tmpli = []
        for i,x in enumerate(intup):
            nolen = len(str(i))
            inds = (5-nolen)*'0'+str(i)
            tmpli.append([x,'SUB'+inds])
        dicli = dict(tmpli)
        dicli2 = dict([[x[1],x[0]] for x in tmpli])
        self.subs_li = deepcopy(tmpli)
        self.var_li = [x[0] for x in tmpli]
        self.subs_dic = deepcopy(dicli)
        self.subs_dic2 = deepcopy(dicli2)

        func = []
        subli = []
        symdic = {}
        for x in dicli.values():
            symdic[x] = sympycore.Symbol(x)
        locals().update(symdic)
#NOTE: don't do this?
# or don't modify _existing_ contents?
# see: http://docs.python.org/library/functions.html#locals
#also locals() is no longer a dictionary in Python 3
#NOTE: just use the dictionary itself
        for x in self.paramdic.keys():
            locals()[x] = sympycore.Symbol(x)
        for x in self.sstate.keys():
            locals()[x] = sympycore.Symbol(x)
        for x in nlsys:
            str_tmp = x[:]
            list1 = self.vreg(patup,x,True,'max')
            list1.reverse()
            for y in list1:
                pos = y[3][0]
                poe = y[3][1]
                vari = y[0]
                str_tmp = str_tmp[:pos] + dicli[vari] +str_tmp[poe:]
            list2 = self.vreg(('{-10,10}|None','iid','{-10,10}'),str_tmp,True,'max')
            if list2:
                list2.reverse()
                for y in list2:
                    pos = y[3][0]
                    poe = y[3][1]
                    vari = y[0]
                    str_tmp = str_tmp[:pos]+'0.0'+str_tmp[poe:]
            # Now substitute out exp and log in terms of sympycore expressions
            elog = re.compile('LOG\(')
            while elog.search(str_tmp):
                ma = elog.search(str_tmp)
                pos = ma.span()[0]
                poe = ma.span()[1]
                str_tmp = str_tmp[:pos]+'sympycore.log('+str_tmp[poe:]
            eexp = re.compile('EXP\(')
            while eexp.search(str_tmp):
                ma = eexp.search(str_tmp)
                pos = ma.span()[0]
                poe = ma.span()[1]
                str_tmp = str_tmp[:pos]+'sympycore.exp('+str_tmp[poe:]
            func.append(eval(str_tmp))
        self.func1 = func
        
        
        # Only exp() when variable needs to be put into ln !
        list1 = [x[1] for x in self.vardic['endo']['var']]
        list1 = list1 + [x[1] for x in self.vardic['con']['var']]
        list1 = list1 + [x[1] for x in self.vardic['exo']['var']]
        list2 = [x for x in self.vardic['endo']['mod']]
        list2 = list2 + [x for x in self.vardic['con']['mod']]
        list2 = list2 + [x for x in self.vardic['exo']['mod']]
        if self.vardic['other']['var']:
            list1 = list1 +  [x[1] for x in self.vardic['other']['var']]
            list2 = list2 +  [x for x in self.vardic['other']['mod']]
        moddic = {}
        for x,y in zip(list1,list2):
            moddic[x] = y
            moddic[x] = y

        lookli = []
        for key in self.vardic.keys():
            lookli = lookli + [[x[0].split('(t)')[0],x[1]] for x in self.vardic[key]['var']]
        lookdic = dict(lookli)
        _mreg = 'SUB\d{5,5}'
        mreg = re.compile(_mreg)
        func2 = []
        for i,x in enumerate(func):
            func2.append(func[i])
            liobj = list(mreg.finditer(str(func2[i])))
            doneli = []
            for ma in reversed(liobj):
                suba = ma.group()
                if 'log' in moddic[lookdic[self.vreg(patup,dicli2[suba],False,'min')[2][1]]]\
                   and suba not in doneli:
                    func2[i] = func2[i].subs(sympycore.Symbol(suba),sympycore.exp(sympycore.Symbol(suba)))
                    doneli.append(suba)
        self.func2 = func2

        # Also take into account log or not for other equations !
        if 'nlsubsys' in dir(self):
            tmpvarli0 = [x[0] for x in self.vardic['other']['var']]
            for i1,line in enumerate(self.func2[-lsubs:]):
                vari0 = tmpvarli0[i1]
                if 'log' in self.vardic['other']['mod'][i1]:
                    pass
                else:
                    vma = self.vreg(patup,vari0,False,'min')
                    vmavar = vma[2][1]+'_bar'
                    self.func2[i1] = sympycore.Symbol(vmavar)*(self.func2[i1])

        # Create list with log of steady states
        evalli = []
        alldic = deepcopy(self.paramdic)
        alldic.update(self.sstate)
        for x in intup:
            vma = self.vreg(patup, x, False, 'min')
            vari = vma[2][1]
            ssvar = float(alldic[vari+'_bar'])
            if 'log' in moddic[lookdic[vari]]:
                evalli.append(np.log(ssvar))
            else:
                evalli.append(ssvar)
        evaldic = {}
        for i,x in enumerate(tmpli):
            evaldic[x[1]] = evalli[i]

        # Make 2D symbolic and numeric Jacobian
        def mkjac(jrows=jrows,jcols=jcols):
            mesg = self._mesg
            rdic = dict([[x,'0'] for x in range(jcols)])
            jdic = dict([[x,rdic.copy()] for x in range(jrows)])
            jdicc = deepcopy(jdic)
            jcols = len(jdic[0])
            jrows = len(jdic)
            carry_over_dic = {}
            numj = mat.zeros((jrows,jcols))
            alldic = {}
            alldic.update(self.paramdic)
            alldic.update(self.sstate)
            alldic.update(evaldic)
            locals().update(alldic)
            for x in range(jrows):
                jdicc[x] = {}
                if mesg:
                    # This will not work in Python 3.0
                    to = "INIT: Computing DSGE model's Jacobian: Equation ("+str(x+1)+"/"+str(jrows)+")..."
                    delete = "\b" * (len(to) + 1)
                    print "{0}{1}".format(delete, to),
                for y in range(jcols):
                    jdic[x][y] = func2[x].diff(symdic[tmpli[y][1]])
                    suba_dic = self.subs_dic2
                    differo_var = str(symdic[tmpli[y][1]])
                    differo_var = suba_dic[differo_var]
                    carry_over_dic[y] = differo_var
                    jdicc[x][differo_var] = str(jdic[x][y])
                    if mreg.search(jdicc[x][differo_var]):
                        for keyo in suba_dic.keys():
                            jdicc[x][differo_var] = jdicc[x][differo_var].replace(keyo,suba_dic[keyo])
                    else:
                        jdicc[x][differo_var] = jdicc[x][differo_var]
                    evalfo = jdic[x][y].evalf()
                    numj[x,y] = eval(str(evalfo))
            # Take out the elements from the variable substitution equations
            lengor = len(self.nlsys_list)
            for keyo in jdicc.keys():
                if keyo > (lengor-1):
                    jdicc.pop(keyo)
            if mesg: print " DONE",
            return numj,jdic,jdicc,carry_over_dic

        # Now make 3D symbolic and numeric Hessian
        def mkhes(jrows=jrows,jcols=jcols,trans_dic=None):
            mesg = self._mesg
            rdic = dict([[x,'0'] for x in range(jcols)])
            rdic1 = dict([[x,rdic.copy()] for x in range(jcols)])
            hdic = dict([[x,deepcopy(rdic1)] for x in range(jrows)])
            hdicc = deepcopy(hdic)
            hcols = len(hdic[0])
            hrows = len(hdic[0])*len(hdic)
            numh = mat.zeros((hrows,hcols))
            jdic = self.jdic
            jdicc = self.jdicc
            alldic = {}
            alldic.update(self.paramdic)
            alldic.update(self.sstate)
            alldic.update(evaldic)
            locals().update(alldic)
            count = 0
            for x in range(jrows):
                hdicc[x] = {}
                if mesg:
                    # This will not work in Python 3.0
                    to = "INIT: Computing DSGE model's Hessian: Equation ("+str(x+1)+"/"+str(jrows)+")..."
                    delete = "\b" * (len(to) + 1)
                    print "{0}{1}".format(delete, to),
                for y in range(jcols):
                    hdicc[x][trans_dic[y]] = {}
                    for z in range(jcols):
                        hdic[x][y][z] = jdic[x][y].diff(symdic[tmpli[z][1]])
                        suba_dic = self.subs_dic2
                        differo_var = str(symdic[tmpli[z][1]])
                        differo_var = suba_dic[differo_var]
                        hdicc[x][trans_dic[y]][differo_var] = str(hdic[x][y][z])
                        if mreg.search(hdicc[x][trans_dic[y]][differo_var]):
                            for keyo in suba_dic.keys():
                                hdicc[x][trans_dic[y]][differo_var] = hdicc[x][trans_dic[y]][differo_var].replace(keyo,suba_dic[keyo])
                        else:
                            hdicc[x][trans_dic[y]][differo_var] = hdicc[x][trans_dic[y]][differo_var]
                        evalfo = hdic[x][y][z].evalf()
                        numh[count,z] = eval(str(evalfo))
                    count = count + 1
            # Take out the elements from the variable substitution equations
            lengor = len(self.nlsys_list)
            for keyo in hdicc.keys():
                if keyo > (lengor-1):
                    hdicc.pop(keyo)
            if mesg: print " DONE",
            return numh,hdic,hdicc

        self.derivatives.numj,self.derivatives.jdic,self.derivatives.jdicc,carry_over_dic = mkjac()
        # To make line between Jacobian's and Hessian's computation
        numj = deepcopy(self.derivatives.numj)
        if mk_hessian:
            self.derivatives.numh,self.derivatives.hdic,self.derivatives.hdicc = mkhes(trans_dic=carry_over_dic)
            numh = deepcopy(self.numh)

        if 'nlsubsys' in dir(self):
            numjs = numj[:-lsubs,:]
            numjl = numj[-lsubs:,:]
            self.derivatives.numj = numjs
            self.derivatives.numjl = numjl
            if mk_hessian:
                numhs = numh[:-lsubs*jcols,:]
                numhl = numh[-lsubs*jcols:,:]
                self.derivatives.numhl = numhl
                self.derivatives.numh = numhs
        else:
            self.derivatives.numj = numj
            if mk_hessian:
                self.derivatives.numh = numh

        self.derivatives.jAA = self.derivatives.numj[:,:int(len(intup)/2)]
        self.derivatives.jBB = -self.derivatives.numj[:,int(len(intup)/2):]
        
        # Get rid of pointless keys in jdicc (and hdicc)
        jdicc_copy = deepcopy(jdicc)
        for keyo in jdicc_copy:
            for keyo2 in jdicc_copy[keyo]:
                if type(keyo2) == type(1): jdicc[keyo].pop(keyo2)
        if mk_hessian:
            hdicc_copy = deepcopy(hdicc)
            if 'hdicc' in dir():
                for keyo in hdicc_copy:
                    for keyo2 in hdicc_copy[keyo]:
                        if type(keyo2) == type(1): hdicc[keyo].pop(keyo2)
                    
        # Build string As and Bs
        sAA = deepcopy(self.derivatives.jAA.astype(str))
        sBB = deepcopy(self.derivatives.jBB.astype(str))
        sAA = np.matrix(sAA,dtype=np.object_)
        sBB = np.matrix(sBB,dtype=np.object_)
        for elem in range(self.derivatives.jAA.shape[0]):
            for i1,elem2 in enumerate(exo_1+endo_1+con_1):
                sAA[elem,i1] = jdicc[elem][elem2]
            for i2,elem2 in enumerate(exo_0+endo_0+con_0):
                if jdicc != '0':
                    sBB[elem,i2] = '-('+jdicc[elem][elem2]+')'
                else:
                    sBB[elem,i2] = jdicc[elem][elem2]     
        self.derivatives.sAA = deepcopy(sAA)
        self.derivatives.sBB = deepcopy(sBB)
        self.derivatives.iA = deepcopy(exo_1+endo_1+con_1)
        self.derivatives.iB = deepcopy(exo_0+endo_0+con_0)
        

    # The parallelized mkjahe version using parallel python

    def mkjahepp(self):
        '''
        A parallelized method using native Python and Sympy in oder
        to calculate the numerical and analytical Jacobian and Hessian of the model.
        This is the parallelized version of method self.mkjahe using the Python pp library.
        
        Parameters
        ----------
        self: object instance
        
        Returns
        -------
        self.numj: attaches numerical model Jacobian to instance
        self.jdic: attaches the analytical model Jacobian to instance
        self.numh: attaches numerical model Hessian to instance
        self.hdic: attaches the analytical 3D Hessian to instance
        self.jAA: attaches numerical AA matrix used in Forkleind solution method
        self.jBB: attaches numerical BB matrix used in Forkleind solution method
        '''
        # Import some configuration options for the DSGE model instance
        ncpus = copy.deepcopy(self._ncpus)
        mk_hessian = copy.deepcopy(self._mk_hessian)
        mesg = copy.deepcopy(self._mesg)
        
        import sympycore
        
        exo_0,exo_1,endo_0,endo_1,con_0,con_1 = self.def_differ_periods()
        
        inlist = exo_1+endo_1+con_1+exo_0+endo_0+con_0
        intup=tuple(inlist)

        patup = ('{-10,10}|None','endo|con|exo|other','{-10,10}')
        jcols = len(intup)
        jrows = len(self.nlsys_list)
        nlsys = copy.deepcopy(self.nlsys_list)

        if 'nlsubsys' in dir(self):
            lsubs = len(self.nlsubsys)
            jrows = jrows + lsubs
            nlsys = nlsys + copy.deepcopy(self.nlsubsys)

        # Create substitution var list and dictionary
        tmpli = []
        for i,x in enumerate(intup):
            nolen = len(str(i))
            inds = (5-nolen)*'0'+str(i)
            tmpli.append([x,'SUB'+inds])
        dicli = dict(tmpli)
        dicli2 = dict([[x[1],x[0]] for x in tmpli])
        self.subs_li = deepcopy(tmpli)
        self.var_li = [x[0] for x in tmpli]
        self.subs_dic = deepcopy(dicli)
        self.subs_dic2 = deepcopy(dicli2)

        func = []
        subli = []
        symdic = {}
        for x in dicli.values():
            symdic[x] = sympycore.Symbol(x)
        locals().update(symdic)
        for x in self.paramdic.keys():
            locals()[x] = sympycore.Symbol(x)
        for x in self.sstate.keys():
            locals()[x] = sympycore.Symbol(x)
        for x in nlsys:
            str_tmp = x[:]
            list1 = self.vreg(patup,x,True,'max')
            list1.reverse()
            for y in list1:
                pos = y[3][0]
                poe = y[3][1]
                vari = y[0]
                str_tmp = str_tmp[:pos] + dicli[vari] +str_tmp[poe:]
            # Take out the IID variables as they don't matter for computation of derivative matrices
            list2 = self.vreg(('{-10,10}|None','iid','{-10,10}'),str_tmp,True,'max')
            if list2:
                list2.reverse()
                for y in list2:
                    pos = y[3][0]
                    poe = y[3][1]
                    vari = y[0]
                    str_tmp = str_tmp[:pos] + '0.0' +str_tmp[poe:]
            # Now substitute out exp and log in terms of sympycore expressions
            elog = re.compile('LOG\(')
            while elog.search(str_tmp):
                ma = elog.search(str_tmp)
                pos = ma.span()[0]
                poe = ma.span()[1]
                str_tmp = str_tmp[:pos]+'sympycore.log('+str_tmp[poe:]
            eexp = re.compile('EXP\(')
            while eexp.search(str_tmp):
                ma = eexp.search(str_tmp)
                pos = ma.span()[0]
                poe = ma.span()[1]
                str_tmp = str_tmp[:pos]+'sympycore.exp('+str_tmp[poe:]
            eexp2 = re.compile('(?<!sympycore\.)E\*\*')
            while eexp2.search(str_tmp):
                ma = eexp2.search(str_tmp)
                pos = ma.span()[0]
                poe = ma.span()[1]
                str_tmp = str_tmp[:pos]+'sympycore.E**'+str_tmp[poe:]
            try:
                func.append(eval(str_tmp))
            except:
                print "ERROR at: "+str_tmp
                sys.exit()
        self.func1 = func

        # Only exp() when variable needs to be put into ln !
        list1 = [x[1] for x in self.vardic['endo']['var']]
        list1 = list1 + [x[1] for x in self.vardic['con']['var']]
        list1 = list1 + [x[1] for x in self.vardic['exo']['var']]
        list2 = [x for x in self.vardic['endo']['mod']]
        list2 = list2 + [x for x in self.vardic['con']['mod']]
        list2 = list2 + [x for x in self.vardic['exo']['mod']]
        if self.vardic['other']['var']:
            list1 = list1 + [x[1] for x in self.vardic['other']['var']]
            list2 = list2 + [x for x in self.vardic['other']['mod']]
        moddic = {}
        for x,y in zip(list1,list2):
            moddic[x] = y
            moddic[x] = y

        lookli = []
        for key in self.vardic.keys():
            lookli = lookli + [[x[0].split('(t)')[0],x[1]] for x in self.vardic[key]['var']]
        lookdic = dict(lookli)
        _mreg = 'SUB\d{5,5}'
        mreg = re.compile(_mreg)
        func2 = []
        for i,x in enumerate(func):
            func2.append(func[i])
            liobj = list(mreg.finditer(str(func2[i])))
            doneli = []
            for ma in reversed(liobj):
                suba = ma.group()
                if 'log' in moddic[lookdic[self.vreg(patup,dicli2[suba],False,'min')[2][1]]]\
                   and suba not in doneli:
                    func2[i] = func2[i].subs(sympycore.Symbol(suba),sympycore.exp(sympycore.Symbol(suba)))
                    doneli.append(suba)
        self.func2 = func2

        # Also take into account log or not for other equations !
        if 'nlsubsys' in dir(self):
            tmpvarli0 = [x[0] for x in self.vardic['other']['var']]
            for i1,line in enumerate(self.func2[-lsubs:]):
                vari0 = tmpvarli0[i1]
                if 'log' in self.vardic['other']['mod'][i1]:
                    pass
                else:
                    vma = self.vreg(patup,vari0,False,'min')
                    vmavar = vma[2][1]+'_bar'
                    self.func2[i1] = sympycore.Symbol(vmavar)*(self.func2[i1])

        # Create list with log of steady states
        evalli = []
        alldic = {}
        alldic.update(self.paramdic)
        alldic.update(self.sstate)
        for x in intup:
            vma = self.vreg(patup, x, False, 'min')
            vari = vma[2][1]
            if 'log' in moddic[lookdic[vari]]:
                evalli.append(np.log(alldic[vari+'_bar']))
            else:
                evalli.append(alldic[vari+'_bar'])
        evaldic = {}
        for i,x in enumerate(tmpli):
            evaldic[x[1]] = evalli[i]

        # Now make the 2D and 3D symbolic and numeric Jacobian and Hessian
        def mkjaheseq(lcount,func2,jcols,symdic,tmpli,paramdic,sstate,evaldic,suba_dic,mk_hessian):
            ### Needed for jdicc and hdicc ###
            _mreg = 'SUB\d{5,5}'
            mreg = re.compile(_mreg)
            ##################################
            
            jdic = dict([[x,'0'] for x in range(jcols)])
            jdicc = copy.deepcopy(jdic)
            carry_over_dic = {}
            numj = numpy.matlib.zeros((1,jcols))
            if mk_hessian:
                rdic = dict([[x,'0'] for x in range(jcols)])
                hdic = dict([[x,rdic.copy()] for x in range(jcols)])
                hdicc = copy.deepcopy(hdic)
                numh = numpy.matlib.zeros((jcols,jcols))

            alldic = {}
            alldic.update(paramdic)
            alldic.update(sstate)
            alldic.update(evaldic)
            locals().update(alldic)
            globals().update(paramdic)
            globals().update(sstate)
            globals().update(evaldic)
            count = 0
            for y in range(jcols):
                jdic[y] = func2[lcount].diff(symdic[tmpli[y][1]])
                
                #### For symbolic jdicc ###
                differo_var = str(symdic[tmpli[y][1]])
                differo_var = suba_dic[differo_var]
                if y not in carry_over_dic.keys(): carry_over_dic[y] = differo_var
                jdicc[differo_var] = str(jdic[y])
                if mreg.search(jdicc[differo_var]):
                    for keyo in suba_dic.keys():
                        jdicc[differo_var] = jdicc[differo_var].replace(keyo,suba_dic[keyo])
                else:
                    jdicc[differo_var] = jdicc[differo_var]
                ###### Done ########
                numj[0,y] = eval(str(jdic[y].evalf()))
                if mk_hessian:
                    for z in range(jcols):
                        hdic[y][z] = jdic[y].diff(symdic[tmpli[z][1]])
                        
                        #### For symbolic hdicc ###
                        differo_var = str(symdic[tmpli[z][1]])
                        differo_var = suba_dic[differo_var]
                        if hdicc.has_key(carry_over_dic[y]):
                            hdicc[carry_over_dic[y]][differo_var] = str(hdic[y][z])
                        else:
                            hdicc[carry_over_dic[y]] = {}
                            hdicc[carry_over_dic[y]][differo_var] = str(hdic[y][z])
                        if mreg.search(hdicc[carry_over_dic[y]][differo_var]):
                            for keyo in suba_dic.keys():
                                hdicc[carry_over_dic[y]][differo_var] = hdicc[carry_over_dic[y]][differo_var].replace(keyo,suba_dic[keyo])
                        else:
                            hdicc[carry_over_dic[y]][differo_var] = hdicc[carry_over_dic[y]][differo_var]
                        ###### Done ########
                        
                        numh[count,z] = eval(str(hdic[y][z].evalf()))
                    count = count + 1
            if mk_hessian:
                return (numj,jdic,jdicc,numh,hdic,hdicc)
            else:
                return (numj,jdic,jdicc)

        inputs = [x for x in xrange(len(self.func2))]
        # Support auto-detection of CPU cores
        if self._ncpus == 'auto':
            if mesg: print "INIT: Parallel execution started with "+str(jobserver.get_ncpus())+ " CPU cores..."
        else:
            if mesg: print "INIT: Parallel execution started with "+str(jobserver.get_ncpus())+ " CPU cores..."

        imports = ('numpy','numpy.matlib','copy','re',)
        
        #job_server.submit(self, func, args=(), depfuncs=(), modules=(), callback=None, callbackargs=(), group='default', globals=None)
        # Submits function to the execution queue
           
        # func - function to be executed
        # args - tuple with arguments of the 'func'
        # depfuncs - tuple with functions which might be called from 'func'
        # modules - tuple with module names to import
        # callback - callback function which will be called with argument
        # list equal to callbackargs+(result,)
        # as soon as calculation is done
        # callbackargs - additional arguments for callback function
        # group - job group, is used when wait(group) is called to wait for
        # jobs in a given group to finish
        # globals - dictionary from which all modules, functions and classes
        # will be imported, for instance: globals=globals()

        # Make sure the jobserver has done his jobs
        jobserver.wait()        
        jobs = [jobserver.submit(mkjaheseq,(inputo,self.func2,jcols,symdic,tmpli,self.paramdic,self.sstate,evaldic,self.subs_dic2,mk_hessian),(),imports) for inputo in inputs]
        if mk_hessian:
            jdic = {}
            jdicc = {}
            hdic = {}
            hdicc = {}
            job_0 = jobs[0]
            numj = job_0()[0]
            jdic[0] = job_0()[1]
            jdicc[0] = job_0()[2]
            numh = job_0()[3]
            hdic[0] = job_0()[4]
            hdicc[0] = job_0()[5]
            for i1,job in enumerate(jobs[1:len(jobs)]):
                numj = mat.vstack((numj,job()[0]))
                jdic[i1+1] = job()[1]
                jdicc[i1+1] = job()[2]
                numh = mat.vstack((numh,job()[3]))
                hdic[i1+1] = job()[4]
                hdicc[i1+1] = job()[5]
            self.derivatives.numj = numj
            self.derivatives.jdic = jdic
            self.derivatives.jdicc = jdicc
            self.derivatives.numh = numh
            self.derivatives.hdic = hdic
            self.derivatives.hdicc = hdicc
        else:
            jdic = {}
            jdicc = {}
            job_0 = jobs[0]
            numj = job_0()[0]
            jdic[0] = job_0()[1]
            jdicc[0] = job_0()[2]
            for i1,job in enumerate(jobs[1:len(jobs)]):
                numj = mat.vstack((numj,job()[0]))
                jdic[i1+1] = job()[1]
                jdicc[i1+1] = job()[2]
            self.derivatives.numj = numj
            self.derivatives.jdic = jdic
            self.derivatives.jdicc = jdicc

        if 'nlsubsys' in dir(self):
            numjs = numj[:-lsubs,:]
            numjl = numj[-lsubs:,:]
            self.derivatives.numj = numjs
            self.derivatives.numjl = numjl

            if mk_hessian:
                numhs = numh[:-lsubs*jcols,:]
                numhl = numh[-lsubs*jcols:,:]
                self.derivatives.numhl = numhl
                self.derivatives.numh = numhs
        else:
            self.numj = numj
            if mk_hessian:
                self.derivatives.numh = numh

        self.derivatives.jAA = self.derivatives.numj[:,:int(len(intup)/2)]
        self.derivatives.jBB = -self.derivatives.numj[:,int(len(intup)/2):]
        
        # Get rid of pointless keys in jdicc (and hdicc)
        jdicc_copy = deepcopy(jdicc)
        for keyo in jdicc_copy:
            for keyo2 in jdicc_copy[keyo]:
                if type(keyo2) == type(1): jdicc[keyo].pop(keyo2)
        if mk_hessian:
            hdicc_copy = deepcopy(hdicc)
            if 'hdicc' in dir():
                for keyo in hdicc_copy:
                    for keyo2 in hdicc_copy[keyo]:
                        if type(keyo2) == type(1): hdicc[keyo].pop(keyo2)
        
        # Build string As and Bs
        sAA = deepcopy(self.derivatives.jAA.astype(str))
        sBB = deepcopy(self.derivatives.jBB.astype(str))
        sAA = np.matrix(sAA,dtype=np.object_)
        sBB = np.matrix(sBB,dtype=np.object_)
        for elem in range(self.derivatives.jAA.shape[0]):
            for i1,elem2 in enumerate(exo_1+endo_1+con_1):
                sAA[elem,i1] = jdicc[elem][elem2]
            for i2,elem2 in enumerate(exo_0+endo_0+con_0):
                if jdicc != '0':
                    sBB[elem,i2] = '-('+jdicc[elem][elem2]+')'
                else:
                    sBB[elem,i2] = jdicc[elem][elem2]     
        self.derivatives.sAA = deepcopy(sAA)
        self.derivatives.sBB = deepcopy(sBB)
        self.derivatives.iA = deepcopy(exo_1+endo_1+con_1)
        self.derivatives.iB = deepcopy(exo_0+endo_0+con_0)        

###NUMERICAL ONLY JACOBIAN AND HESSIAN METHODS - FOR UPDATING###############
    # The unparallelized Jacobian, Hessian computation method
    def mkjahen(self):
        '''
        An unparallelized method using native Python and Sympycore in oder
        to calculate the numerical and analytical Jacobian and Hessian of the model.
        This is the corresponding method to "self.mkjahe()" but is only called when
        some parameters were changed and the model needs dynamic updating. In particular
        in this case the expensive computation of the analytical expressions can be
        avoided, as only parameters affecting the steady state may have been altered.
        Uses existing self.jdic for Jacobian and self.hdic for Hessian of model.
        
        Parameters
        ----------
        self: object instance
        
        Returns
        -------
        self.numj: attaches numerical model Jacobian to instance
        self.numh: attaches numerical model Hessian to instance
        self.jAA: attaches numerical AA matrix used in Forkleind solution method
        self.jBB: attaches numerical BB matrix used in Forkleind solution method
        '''
        mk_hessian = self._mk_hessian
        
        exo_0,exo_1,endo_0,endo_1,con_0,con_1 = self.def_differ_periods()
        inlist = exo_1+endo_1+con_1+exo_0+endo_0+con_0
        intup=tuple(inlist)

        patup = ('{-10,10}|None','endo|con|exo|other','{-10,10}')
        jcols = len(intup)
        jrows = len(self.nlsys_list)
        nlsys = deepcopy(self.nlsys_list)

        if 'nlsubsys' in dir(self):
            lsubs = len(self.nlsubsys)
            jrows = jrows + lsubs
            nlsys = nlsys + deepcopy(self.nlsubsys)

        # Create substitution var list and dictionary
        tmpli = []
        for i,x in enumerate(intup):
            nolen = len(str(i))
            inds = (5-nolen)*'0'+str(i)
            tmpli.append([x,'SUB'+inds])
            dicli = dict(tmpli)
            dicli2 = dict([[x[1],x[0]] for x in tmpli])
            
        lookli = []
        for key in self.vardic.keys():
            lookli = lookli + [[x[0].split('(t)')[0],x[1]] for x in self.vardic[key]['var']]
        lookdic = dict(lookli)

        # Only exp() when variable needs to be put into ln !
        list1 = [x[1] for x in self.vardic['endo']['var']]
        list1 = list1 + [x[1] for x in self.vardic['con']['var']]
        list1 = list1 + [x[1] for x in self.vardic['exo']['var']]
        list2 = [x for x in self.vardic['endo']['mod']]
        list2 = list2 + [x for x in self.vardic['con']['mod']]
        list2 = list2 + [x for x in self.vardic['exo']['mod']]
        if self.vardic['other']['var']:
            list1 = list1 +  [x[1] for x in self.vardic['other']['var']]
            list2 = list2 +  [x for x in self.vardic['other']['mod']]
        moddic = {}
        for x,y in zip(list1,list2):
            moddic[x] = y
            moddic[x] = y

        # Create list with log of steady states
        evalli = []
        alldic = {}
        alldic.update(self.paramdic)
        alldic.update(self.sstate)
        for x in intup:
            vma = self.vreg(patup, x, False, 'min')
            vari = vma[2][1]
            if 'log' in moddic[lookdic[vari]]:
                evalli.append(np.log(alldic[vari+'_bar']))
            else:
                evalli.append(alldic[vari+'_bar'])
        evaldic = {}
        for i,x in enumerate(tmpli):
            evaldic[x[1]] = evalli[i]

        # Make 2D symbolic and numeric Jacobian
        def mkjac(jrows=jrows,jcols=jcols):
            jdic = self.jdic
            numj = mat.zeros((jrows,jcols))
            alldic = {}
            alldic.update(self.paramdic)
            alldic.update(self.sstate)
            alldic.update(evaldic)
            locals().update(alldic)
            for x in range(jrows):
                for y in range(jcols):
                    evalfo = jdic[x][y].evalf()
                    if 'exp(' not in str(evalfo) and 'log(' not in str(evalfo):
                        numj[x,y] = eval(str(evalfo))
                    elif 'exp(' in str(evalfo):
                        numj[x,y] = eval(str(evalfo).replace('exp(','np.exp('))
                    elif 'log(' in str(evalfo):
                        numj[x,y] = eval(str(evalfo).replace('log(','np.log('))
            return numj

        # Now make 3D symbolic and numeric Hessian
        def mkhes(jrows=jrows,jcols=jcols):
            hdic = self.hdic
            hrows = jrows*jcols
            numh = mat.zeros((hrows,jcols))
            alldic = {}
            alldic.update(self.paramdic)
            alldic.update(self.sstate)
            alldic.update(evaldic)
            locals().update(alldic)
            count = 0
            for x in range(jrows):
                for y in range(jcols):
                    for z in range(jcols):
                        evalfo = hdic[x][y][z].evalf()
                        if 'exp(' not in str(evalfo) and 'log(' not in str(evalfo):
                            numh[count,z] = eval(str(evalfo))
                        elif 'exp(' in str(evalfo):
                            numh[count,z] = eval(str(evalfo).replace('exp(','np.exp('))
                        elif 'log(' in str(evalfo):
                            numh[count,z] = eval(str(evalfo).replace('log(','np.log('))
                    count = count + 1
            return numh

        self.derivatives.numj = mkjac()
        numj = deepcopy(self.numj)
        if mk_hessian:
            self.derivatives.numh = mkhes()
            numh = deepcopy(self.numh)

        if 'nlsubsys' in dir(self):
            numjs = numj[:-lsubs,:]
            numjl = numj[-lsubs:,:]
            self.derivatives.numj = numjs
            self.derivatives.numjl = numjl
            if mk_hessian:
                numhs = numh[:-lsubs*jcols,:]
                numhl = numh[-lsubs*jcols:,:]
                self.derivatives.numhl = numhl
                self.derivatives.numh = numhs
        else:
            self.derivatives.numj = numj
            if mk_hessian:
                self.derivatives.numh = numh

        self.derivatives.jAA = self.derivatives.numj[:,:int(len(intup)/2)]
        self.derivatives.jBB = -self.derivatives.numj[:,int(len(intup)/2):]

    def mkjaheppn(self):
        '''
        A parallelized method using native Python and Sympycore in oder
        to calculate the numerical and analytical Jacobian and Hessian of the model.
        This is the corresponding method to "self.mkjahe()" but is only called when
        some parameters were changed and the model needs dynamic updating. In particular
        in this case the expensive computation of the analytical expressions can be
        avoided, as only parameters affecting the steady state may have been altered.
        Uses existing self.jdic for Jacobian and self.hdic for Hessian of model.
        
        Parameters
        ----------
        self: object instance
        
        Returns
        -------
        self.numj: attaches numerical model Jacobian to instance
        self.numh: attaches numerical model Hessian to instance
        self.jAA: attaches numerical AA matrix used in Forkleind solution method
        self.jBB: attaches numerical BB matrix used in Forkleind solution method
        '''  

        mk_hessian = self._mk_hessian
        # import local sympycore
        import sympycore
        
        exo_0,exo_1,endo_0,endo_1,con_0,con_1 = self.def_differ_periods()
        inlist = exo_1+endo_1+con_1+exo_0+endo_0+con_0
        intup=tuple(inlist)

        patup = ('{-10,10}|None','endo|con|exo|other','{-10,10}')
        jcols = len(intup)
        jrows = len(self.nlsys_list)
        nlsys = deepcopy(self.nlsys_list)

        if 'nlsubsys' in dir(self):
            lsubs = len(self.nlsubsys)
            jrows = jrows + lsubs
            nlsys = nlsys + deepcopy(self.nlsubsys)

        # Create substitution var list and dictionary
        tmpli = []
        for i,x in enumerate(intup):
            nolen = len(str(i))
            inds = (5-nolen)*'0'+str(i)
            tmpli.append([x,'SUB'+inds])
            dicli = dict(tmpli)
            dicli2 = dict([[x[1],x[0]] for x in tmpli])

        lookli = []
        for key in self.vardic.keys():
            lookli = lookli + [[x[0].split('(t)')[0],x[1]] for x in self.vardic[key]['var']]
        lookdic = dict(lookli)

        # Only exp() when variable needs to be put into ln !
        list1 = [x[1] for x in self.vardic['endo']['var']]
        list1 = list1 + [x[1] for x in self.vardic['con']['var']]
        list1 = list1 + [x[1] for x in self.vardic['exo']['var']]
        list2 = [x for x in self.vardic['endo']['mod']]
        list2 = list2 + [x for x in self.vardic['con']['mod']]
        list2 = list2 + [x for x in self.vardic['exo']['mod']]
        if self.vardic['other']['var']:
            list1 = list1 +  [x[1] for x in self.vardic['other']['var']]
            list2 = list2 +  [x for x in self.vardic['other']['mod']]
        moddic = {}
        for x,y in zip(list1,list2):
            moddic[x] = y
            moddic[x] = y

        # Create list with log of steady states
        evalli = []
        alldic = {}
        alldic.update(self.paramdic)
        alldic.update(self.sstate)
        for x in intup:
            vma = self.vreg(patup, x, False, 'min')
            vari = vma[2][1]
            if 'log' in moddic[lookdic[vari]]:
                evalli.append(np.log(alldic[vari+'_bar']))
            else:
                evalli.append(alldic[vari+'_bar'])
        evaldic = {}
        for i,x in enumerate(tmpli):
            evaldic[x[1]] = evalli[i]

        # Now make the 2D and 3D symbolic and numeric Jacobian and Hessian
        def mkjaheseq(lcount,func2,jcols,tmpli,paramdic,sstate,evaldic,mk_hessian,jdic,hdic):
            numj = numpy.matlib.zeros((1,jcols))
            if mk_hessian:
                numh = numpy.matlib.zeros((jcols,jcols))

            alldic = {}
            alldic.update(paramdic)
            alldic.update(sstate)
            alldic.update(evaldic)
            locals().update(alldic)
            globals().update(paramdic)
            globals().update(sstate)
            globals().update(evaldic)
            count = 0
            for y in range(jcols):
                evalfo = jdic[lcount][y].evalf()
                if 'exp(' not in str(evalfo) and 'log(' not in str(evalfo):
                    numj[0,y] = eval(str(evalfo))
                elif 'exp(' in str(evalfo):
                    numj[0,y] = eval(str(evalfo).replace('exp(','np.exp('))
                elif 'log(' in str(evalfo):
                    numj[0,y] = eval(str(evalfo).replace('log(','np.log('))
                if mk_hessian:
                    for z in range(jcols):
                        evalfo2 = hdic[lcount][y][z].evalf()
                        if 'exp(' not in str(evalfo2) and 'log(' not in str(evalfo2):
                            numh[count,z] = eval(str(evalfo2))
                        elif 'exp(' in str(evalfo2):
                            numh[count,z] = eval(str(evalfo2).replace('exp(','np.exp('))
                        elif 'log(' in str(evalfo2):
                            numh[count,z] = eval(str(evalfo2).replace('log(','np.log('))
                    count = count + 1
            if mk_hessian:
                return (numj,numh)
            else:
                return numj


        inputs = [x for x in xrange(len(self.func2))]
        if self._ncpus == 'auto':
            if self._mesg: print "INIT: Parallel execution started with "+str(jobserver.get_ncpus())+ " CPU cores..."
        else:
            if self._mesg: print "INIT: Parallel execution started with "+str(jobserver.get_ncpus())+ " CPU cores..."
        imports = ('numpy','copy','numpy.matlib',)
        # Make sure the jobserver has done his jobs
        jobserver.wait()        
        jobs = [jobserver.submit(mkjaheseq,(inputo,self.func2,jcols,tmpli,self.paramdic,self.sstate,evaldic,mk_hessian,self.jdic,self.hdic),(),imports) for inputo in inputs]
        if mk_hessian:
            job_0 = jobs[0]
            numj = job_0()[0]
            numh = job_0()[1]
            for i1,job in enumerate(jobs[1:len(jobs)]):
                numj = mat.vstack((numj,job()[0]))
                numh = mat.vstack((numh,job()[1]))
            self.derivatives.numj = numj
            self.derivatives.numh = numh
        else:
            job_0 = jobs[0]
            numj = job_0()
            for i1,job in enumerate(jobs[1:len(jobs)]):
                numj = mat.vstack((numj,job()))
            self.derivatives.numj = numj

        if 'nlsubsys' in dir(self):
            numjs = numj[:-lsubs,:]
            numjl = numj[-lsubs:,:]
            self.derivatives.numj = numjs
            self.derivatives.numjl = numjl
            
            if mk_hessian:
                numhs = numh[:-lsubs*jcols,:]
                numhl = numh[-lsubs*jcols:,:]
                self.derivatives.numhl = numhl
                self.derivatives.numh = numhs
        else:
            self.derivatives.numj = numj
            if mk_hessian:
                self.derivatives.numh = numh

        self.derivatives.jAA = self.derivatives.numj[:,:int(len(intup)/2)]
        self.derivatives.jBB = -self.derivatives.numj[:,int(len(intup)/2):]       
##########################################################  
    # Method updating model IF model file has been changed externally !
    def updf(self):
        '''
        Method useful for completely re-initializing model based on a different or else
        altered external modfile. It is also perceivable to load the existing modfile at
        runtime, change it in a program, save it as a temporary file and then re-load it
        using this method.
        '''
        self.txtpars.__init__(self.modfile)


