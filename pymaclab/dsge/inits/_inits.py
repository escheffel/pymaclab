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
import pp

# Import of refactored steady state solution methods
from pymaclab.dsge.solvers.steadystate import SSsolvers, ManualSteadyState, Fsolve

# Import of refactored basic pymaclab parsing method
from pymaclab.dsge.parsers._modparser import parse_mod

# Import of refactored dsge_parser functions to which DSGE instance will be passed
from pymaclab.dsge.parsers._dsgeparser import populate_model_stage_one,populate_model_stage_one_a,\
     populate_model_stage_one_b,populate_model_stage_one_bb,populate_model_stage_two

# Import of refactored updater Classes/methods
from pymaclab.dsge.updaters.one_off import Updaters, dicwrap, dicwrap_deep, listwrap, matwrap
from pymaclab.dsge.updaters.queued import Updaters_Queued, dicwrap_queued, dicwrap_deep_queued,\
     listwrap_queued, matwrap_queued, Process_Queue, queue

# Import of Translator Class for translating between different model formats
from pymaclab.dsge.translators.translators import Translators

# Import refactored DSGE model initialisor class and the dynarepp flag
import pymaclab.dsge.inits._set_flags
from pymaclab.dsge.inits._set_flags import dynare_flag

# Import refactored additional helpful objects
from pymaclab.dsge.helpers._paths import *
from pymaclab.dsge.helpers._greek_alphab import greek_alph

# Instantiate a global class jobserver accessible from all instances
global ppservers, jobserver
ppservers = ()
jobserver = pp.Server(ppservers=ppservers)


class Inits(object):
    '''
    This is the factored-out macrolab DSGEmodel Initialisation class. This is used to be part of the DSGEmodel
    class itself in the form of methods, but has been factored out into a separate class which takes a newly
    instantiated DSGEmodel instance as and argument. It then possess several init methods which "work on the
    instance" and produce desired attributes to different depths of solution levels.
    
    .. note::
    
       Notice that the various init method, i.e. init1, init1a, etc. are not called here, but they are called
       externally in the pymaclab package when instantiating a new DSGE model using pymaclab.newMOD().
    
    :param other:             The DSGEmodel instance to be passed to the Inits instance
    :type other:              DSGEmodel

    '''
    def __init__(self,other=None):
        self._other = other

    def init1(self,no_wrap=False):
        '''
        The init1 method. Model population proceeds from the __init__function here. In particular the data gets read in
        (not implemented at the moment) and the model parsing begins.
        
        .. note::
           The other.vardic gets created and the manual as well as
           the numerical steady state sections gets parsed and attached to the DSGE model instance. So the most import fields created
           here using function populate_model_stage_one() are:
           
             * other.vardic - variable names dictionary with transform and filtering info
             * other.mod_name - short name of the DSGE model from mod file
             * other.mod_desc - longer model description from mod file
             * other.paramdic - dic of defined parameters with their numerical values
             * other.manss_sys - list of equations from the closed form steady state section
             * other.ssys_list - list of equations from the numerical steady state section
           
           Also the updaters and updaters_queued branches are opened here and the other.vardic gets wrapped for dynamic updating
           behaviour.
        
        :param other:     The DSGE model instance itother.
        :type other:      dsge_inst
        
        '''
        other = self._other
        initlev = other._initlev
        ncpus = other._ncpus
        mk_hessian = other._mk_hessian
        mesg = other._mesg
        # Attach the data from database
        if other.dbase != None:
            other.getdata(dbase=other.dbase)        
        if mesg: print "INIT: Instantiating DSGE model with INITLEV="+str(initlev)+" and NCPUS="+str(ncpus)+"..."
        if mesg: print "INIT: Attaching model properties to DSGE model instance..."

        # Create None tester regular expression
        #        _nreg = '^\s*None\s*$'
        #        nreg = re.compile(_nreg)
        
        if mesg: print "INIT: Parsing model file into DSGE model instance..."
        txtpars = parse_mod(other.modfile)
        other.txtpars = txtpars  # do we need txtpars attached for anything else?
        secs = txtpars.secs # do we need txtpars attached for anything else?
        if mesg: print "INIT: Extraction of info into DSGE model instance Stage [1]..."
        # Initial population method of model, does NOT need steady states
        if txtpars.fmt == 'pymaclab':
            other = populate_model_stage_one(other, secs)
        elif txtpars.fmt == 'dynarepp':
            pml_modfile = dynarepp_to_pml.translate(secs=secs)

        if not no_wrap:
            # Open updaters path
            other.updaters = Updaters()        
            # Open the updaters_queued path
            other.updaters_queued = Updaters_Queued()
            # Add the empty queue
            other.updaters_queued.queue = queue
    
            # Wrap the vardic
            other.updaters.vardic = dicwrap_deep(other,'self.vardic',initlev)
            other.updaters_queued.vardic = dicwrap_deep_queued(other,'self.vardic',initlev)
################################### BASIC READING-IN OF MODEL FILE INFO DONE ###############################        
    def init1a(self):
        '''
        The init1a method. Model population proceeds from the init1 method here. The only field which get created here is the raw
        (i.e. unsubstituted) substitution dictionary.
        
        .. note::
            Field which are created here using the function populate_model_stage_one_a() which in turn calls mk_subs_dic():
           
             * other.nlsubs_raw1 - a list of the @items and their replacements
             * other.nlsubsdic - the above just expressed as a keyed list
        
        :param other:     The DSGE model instance itother.
        :type other:      dsge_inst

        '''
        other = self._other
        mesg = other._mesg
        secs = other.txtpars.secs
        # This is an additional populator which creates subsitution dictionary
        # and uses this to already get rid of @s in the manuall sstate list
        if mesg: print "INIT: Extraction of info into DSGE model instance Stage [2]..."
        # This function only creates the raw substitution dictionary and list from modfile
        other = populate_model_stage_one_a(other,secs)
        
    def init1b(self,no_wrap=False):
        '''
        The init1b method. Model population proceeds from the init1a method here.

        .. note::
        
           The only thing which gets done here purposefully *after* calling init1a() is to wrap these fields:
           
             * other.nlsubsdic - the above just expressed as a keyed list
             * other.paramdic - dic of defined parameters with their numerical values
             
           in order to give them dynamic updating behaviour. No more is done in this init method call.
        
        :param other:     The DSGE model instance itother.
        :type other:      dsge_inst

        '''
        other = self._other
        initlev = other._initlev
        secs = other.txtpars.secs
        other = populate_model_stage_one_b(other,secs)
        
        if not no_wrap:
            # Wrap the nlsubsdic
            if 'nlsubsdic' in dir(other): other.updaters_queued.nlsubsdic = dicwrap_queued(other,'self.nlsubsdic',initlev)
            # Wrap the paramdic
            other.updaters_queued.paramdic = dicwrap_queued(other,'self.paramdic',initlev)
            
    
            # Wrap the nlsubsdic
            if 'nlsubsdic' in dir(other): other.updaters.nlsubsdic = dicwrap(other,'self.nlsubsdic',initlev)
            # Wrap the paramdic
            other.updaters.paramdic = dicwrap(other,'self.paramdic',initlev)        
        
    def init1c(self,no_wrap=False):
        '''
        The init1c method. Model population proceeds from the init1b method here. We call populate_model_stage_one_bb() which does quite
        a bit of substitution/replacement of the @-prefixed variables. It does this in the numerical and closed form steady state
        calculation sections, as well as in the actual FOCs themselves.

        .. note::
        
           *After* that we can wrap various fields for updating which are:
           
             * other.foceqs - list of firs-order conditions with @ replacements done
             * other.manss_sys - the list of equations from closed form SS section
             * other.syss_list - the list of equations from numerical SS section
        
           Notice also that here we replace or generate fields in case the FOCs are supposed to be used
           directly in SS calculation because the "use_focs" parameter was not passed empty.
        
        :param other:     The DSGE model instance itother.
        :type other:      dsge_inst

        '''
        other = self._other
        initlev = other._initlev
        secs = other.txtpars.secs
        mesg = other._mesg
        if mesg: print "INIT: Substituting out @ variables in steady state sections..."
        other = populate_model_stage_one_bb(other,secs)
        
        if not no_wrap:
            if 'foceqs' in dir(other):
                # Wrap foceqs
                other.updaters.foceqs = listwrap(other,'self.foceqs',initlev)
                # Wrap foceqs
                other.updaters_queued.foceqs = listwrap(other,'self.foceqs',initlev)
        
        # Allow for numerical SS to be calculated using only FOCs
        if other._use_focs and other._ssidic != None:
            list_tmp = []
            for elem in other._use_focs:
                list_tmp.append(other.foceqss[elem])
            other.ssys_list = deepcopy(list_tmp)
            # Check if this has not already been produced in populate_model_stage_one_bb, chances are it has
            if 'ssidic' not in dir(other): other.ssidic = copy.deepcopy(other._ssidic)
        
        if not no_wrap:   
            # Wrap manss_sys
            if 'manss_sys' in dir(other):
                other.updaters.manss_sys = listwrap(other,'self.manss_sys',initlev)
                other.updaters_queued.manss_sys = listwrap_queued(other,'self.manss_sys',initlev)
            # Wrap ssys_list
            if 'ssys_list' in dir(other):
                other.updaters.ssys_list = listwrap(other,'self.ssys_list',initlev)
                other.updaters_queued.ssys_list = listwrap_queued(other,'self.ssys_list',initlev)


    def init2(self):
        '''
        The init2 method. Model population proceeds from the init1c method here. In this initialisation method
        the only thing which is being done is to open up the sssolvers branch and pass down required objects to
        the manuals closed from solver and the numerical root-finding solver depending on whether information for
        this has been included in the DSGE model file. *No* attempt is made at solving for the steady state, the
        respective solvers are only being *prepared*.
        
        :param other:     The DSGE model instance itother.
        :type other:      dsge_inst
        
        '''
        other = self._other
        mesg = other._mesg
        secs = other.txtpars.secs
        initlev = other._initlev
        if mesg: print "INIT: Preparing DSGE model instance for steady state solution..."

        # Attach the steady state class branch, and add Fsolve if required but do no more !
        other.sssolvers = SSsolvers()

        # check if the Steady-State Non-Linear system .mod section has an entry
        if all([False if 'None' in x else True for x in secs['manualss'][0]]) or other._use_focs:
            intup = (other.ssys_list,other.ssidic,other.paramdic)
            other.sssolvers.fsolve = Fsolve(intup,other=other)
        # check if the Steady-State Non-Linear system closed form has an entry
        #  we do this here with ELIF because we only want to set this up for solving if manualss does not exist
        elif all([False if 'None' in x else True for x in secs['closedformss'][0]]):
            alldic = {}
            alldic.update(other.paramdic)
            intup = (other.manss_sys,alldic)
            other.sssolvers.manss = ManualSteadyState(intup,other=other) 
##################################### FULL PARSING AND UPDATER WRAPPING DONE #############################
    def init3(self):
        '''
        Initialisation method init3 is quite complex and goes through a number of logical tests in order to determine
        how to solve for the steady state of the model.
        
        .. note::
        
           There are 7 different ways a DSGE model can obtain its steady state solution depending on what information has
           been provided:
           
             1) The steady state values dictionary has been passed as argument, then init3() will NEVER be called
             2) Information has been provided using the "use_focs" parameter to use FOCs directly, externally passed using use_focs
             3) Information has been provided using the "use_focs" parameter to use FOCs directly, but inside model file
             4) Information has only been provided in the numerical SS section
             5) Information has only been provided in the closed form SS section
             6) Both CF-SS and NUM-SS info are present and NUM-SS is subset if CF-SS
             7) Both CF-SS and NUM-SS info are present and CF is residual
             
           These options are better explained in the documentation to PyMacLab in the steady state solver section.
           
        :param other:     The DSGE model instance itother.
        :type other:      dsge_inst
        
        '''
        other = self._other
        txtpars = other.txtpars
        secs = txtpars.secs
        initlev = other._initlev
################################## STEADY STATE CALCULATIONS !!! #######################################
        if other._mesg: print "INIT: Attempting to find DSGE model's steady state automatically..."
        # ONLY NOW try to solve !
        ##### OPTION 1: There is only information externally provided and we are using FOCs
        if other._use_focs and other._ssidic != None:
            if '_internal_focs_used' not in dir(other):
                if other._mesg: print "SS: Using FOCs and EXTERNALLY supplied information...attempting to solve SS..."
            else:
                ##### OPTION 1b: There is only information internally provided and we are using FOCs
                if other._mesg: print "SS: Using FOCs and INTERNALLY supplied information...attempting to solve SS..."
                del other._internal_focs_used
            other.sssolvers.fsolve.solve()
            if other.sssolvers.fsolve.ier == 1:
                other.sstate = deepcopy(other.sssolvers.fsolve.fsout)
                other.switches['ss_suc'] = ['1','1']
                if other._mesg: print "INIT: Steady State of DSGE model found (SUCCESS)..."
            else:
                other.switches['ss_suc'] = ['1','0']
                if other._mesg: print "INIT: Steady State of DSGE model not found (FAILURE)..."            
            return
            
        ##### OPTION 2: There is only information provided in the numerical section NOT in closed-form
        if not all([False if 'None' in x else True for x in secs['closedformss'][0]]) and\
           all([False if 'None' in x else True for x in secs['manualss'][0]]) and not other._use_focs and other._ssidic == None:
            if other._mesg: print "SS: ONLY numerical steady state information supplied...attempting to solve SS..."
            other.sssolvers.fsolve.solve()

        ##### OPTION 3: There is only information on closed-form steady state, BUT NO info on numerical steady state
        # Solve using purely closed form solution if no other info on model is available
        if all([False if 'None' in x else True for x in secs['closedformss'][0]]) and\
           not all([False if 'None' in x else True for x in secs['manualss'][0]]):
            if other._mesg: print "SS: ONLY CF-SS information supplied...attempting to solve SS..."
            alldic = {}
            alldic.update(other.paramdic)
            intup = (other.manss_sys,alldic)
            other.sssolvers.manss = ManualSteadyState(intup)
            other.sssolvers.manss.solve()
            other.sstate = {}
            other.sstate.update(other.sssolvers.manss.sstate)
            # Save the solution of the closed-form calculations in the template_paramdic for dynare
            other.template_paramdic['ssidic'] = {}
            other.template_paramdic['ssidic'].update(other.sstate)
            # Also save as ssili for consistency's sake
            other.template_paramdic['ssili'] = []
            for itemo in other.template_paramdic['ssidic'].items():
                other.template_paramdic['ssili'].append([itemo[0],str(itemo[1])])
        ##### OPTION 4: There is information on closed-form AND on numerical steady state
        # Check if the numerical and closed form sections have entries
        if all([False if 'None' in x else True for x in secs['manualss'][0]]) and\
           all([False if 'None' in x else True for x in secs['closedformss'][0]]):
            # Create unordered Set of closed from solution variables
            manss_set = set()
            for elem in other.manss_sys:
                manss_set.add(elem.split('=')[0].lstrip().rstrip())
            # If ssidic is not empty we need to make sure it perfectly overlaps with manss_set in order to replace ssidic
            if other.ssidic != {}:
                numss_set = set()
                for keyo in other.ssidic.keys():
                    numss_set.add(keyo)
                ##### OPTION 4a: If there is an ssidic and its keys are subset of manss_set, the use as suggestion for new ssi_dic
                if numss_set.issubset(manss_set):
                    if other._mesg: print "SS: CF-SS and NUM-SS (overlapping) information information supplied...attempting to solve SS..."
                    alldic = {}
                    alldic.update(other.paramdic)
                    intup = (other.manss_sys,alldic)
                    other.sssolvers.manss = ManualSteadyState(intup)
                    other.sssolvers.manss.solve()
                    for keyo in other.ssidic.keys():
                        other.ssidic[keyo] = other.sssolvers.manss.sstate[keyo]
            ##### OPTION 4b: ssidic is empty, so we have to assumed that the variables in closed form are suggestions for ssidic
            # If it is empty, then just compute the closed form SS and pass to ssidic as starting value
            elif other.ssidic == {}:
                if other._mesg: print "SS: CF-SS and NUM-SS (empty ssdic) information information supplied...attempting to solve SS..."
                alldic = {}
                alldic.update(other.paramdic)
                intup = (other.manss_sys,alldic)
                other.sssolvers.manss = ManualSteadyState(intup)
                other.sssolvers.manss.solve()
                other.ssidic.update(other.sssolvers.manss.sstate)
                # Test at least if the number of ssidic vars equals number of equations
                if len(other.ssidic.keys()) != len(other.ssys_list):
                    print "Error: Number of variables in initial starting values dictionary != number of equations"
                    sys.exit()
        ######## Finally start the numerical root finder with old or new ssidic from above
        if all([False if 'None' in x else True for x in secs['manualss'][0]]) and\
           not all([False if 'None' in x else True for x in secs['closedformss'][0]]) and not other._use_focs:            
            other.sssolvers.fsolve.solve()
        elif all([False if 'None' in x else True for x in secs['manualss'][0]]) and\
           all([False if 'None' in x else True for x in secs['closedformss'][0]]) and not other._use_focs:           
            other.sssolvers.fsolve.solve()

        if all([False if 'None' in x else True for x in secs['manualss'][0]]):
            if other.sssolvers.fsolve.ier == 1:
                other.sstate = other.sssolvers.fsolve.fsout
                other.numssdic = other.sssolvers.fsolve.fsout
                # Attach solutions to intial variable dictionaries, for further analysis
                other.ssidic_modfile = deepcopy(other.ssidic)
                # Update old ssidic with found solutions
                other.ssidic = deepcopy(other.sssolvers.fsolve.fsout)
                other.sssolvers.fsolve.ssi = other.ssidic
                other.switches['ss_suc'] = ['1','1']
                if other._mesg: print "INIT: Steady State of DSGE model found (SUCCESS)..."
            else:
                other.switches['ss_suc'] = ['1','0']
                if other._mesg: print "INIT: Steady State of DSGE model not found (FAILURE)..."

        ########## Here we are trying to merge numerical SS solver's results with result closed-form calculations, if required
        if all([False if 'None' in x else True for x in secs['manualss'][0]]) and\
           all([False if 'None' in x else True for x in secs['closedformss'][0]]):
            if other._mesg: print "INIT: Merging numerical with closed form steady state if needed..."
        # Check for existence of closedform AND numerical steady state
        # We need to stop the model instantiation IFF numerical solver was attempted but failed AND closed form solver depends on it.
        if all([False if 'None' in x else True for x in secs['closedformss'][0]]) and\
           all([False if 'None' in x else True for x in secs['manualss'][0]]) and other.ssidic != {}:
            if other.switches['ss_suc'] == ['1','0']:
                print "ERROR: You probably want to use numerical steady state solution to solve for RESIDUAL closed form steady states."
                print "However, the numerical steady state solver FAILED to find a root, so I am stopping model instantiation here."
                sys.exit()
            ##### OPTION 5: We have both numerical and (residual) closed form information
            # Check if a numerical SS solution has been attempted and succeeded, then take solutions in here for closed form.
            elif other.switches['ss_suc'] == ['1','1'] and not numss_set.issubset(manss_set):
                if other._mesg: print "SS: CF-SS (residual) and NUM-SS information information supplied...attempting to solve SS..."
                alldic = {}
                alldic.update(other.sstate)
                alldic.update(other.paramdic)
                intup = (other.manss_sys,alldic)
                other.sssolvers.manss = ManualSteadyState(intup)
                other.sssolvers.manss.solve()
                # Here merging takes place
                other.sstate.update(other.sssolvers.manss.sstate)

        # Double check if no steady state values are negative, as LOGS may have to be taken.
        if 'sstate' in dir(other):
            for keyo in other.sstate.keys():
                if '_bar' in keyo and float(other.sstate[keyo]) < 0.0:
                    print "WARNING: Steady state value "+keyo+ " is NEGATIVE!"
                    print "This is very likely going to either error out or produce strange results"
                    print "Re-check your model declarations carefully!"
        # Re-attach the solutions back to the dynare template generator dictionary
        for keyo in other.template_paramdic['ssidic'].keys():
            if keyo in other.sstate.keys(): other.template_paramdic['ssidic'][keyo] = deepcopy(other.sstate[keyo])
            other.template_paramdic['ssili'][[x[0] for x in other.template_paramdic['ssili']].index(keyo)][1] = str(other.sstate[keyo])
##################################### STEADY STATE CALCULATION SECTION DONE ##############################
    def init4(self,no_wrap=False):
        '''
        This model instance sub-initializor only calls the section which use the computed steady state
        in order to compute derivatives and open dynamic solver branches on the instance. But Jacobian and Hessian
        are *not* computed here, this is postponed to the next init level.
        
        .. note::
        
           The following last field is wrapped for dynamic execution:
           
             * other.sigma - the variance-covariance matrix of the iid shocks
             
           Notice that after wrapping this last field the process_queue class is instantiated at last,
           because it needs to have access to *all* of the wrapped fields. Also in this method, the function
           populate_model_stage_two() is called which prepares the nonlinear FOCs for derivative-taking.
           
        :param other:     The DSGE model instance itother.
        :type other:      dsge_inst

        '''
        other = self._other
        txtpars = other.txtpars
        secs = txtpars.secs
        initlev = other._initlev
        ncpus = other._ncpus
        mk_hessian = other._mk_hessian
        mesg = other._mesg
        if mesg: print "INIT: Preparing DSGE model instance for computation of Jacobian and Hessian..."
        # Now populate more with stuff that needs steady state
        other = populate_model_stage_two(other, secs)
        # We can also already compute a dictionary holding all the vars w.r.t which derivatives will be taken
        # This will attach a dictionary to the model instance called self.differvardic
        other.def_differ_periods()
        # Open the translators branch on the DSGE model instance
        # This should be done here, so that the passed "other" model instance has the translator branch
        other.translators = Translators(other=other)         

        if not no_wrap:
            # Need to wrap variance covariance matrix here
            other.updaters.sigma = matwrap(other,'self.sigma',initlev)
            # Need to wrap variance covariance matrix here
            other.updaters_queued.sigma = matwrap_queued(other,'self.sigma',initlev)

            ####### All queued updaters initialized, no add processing instance
            # Add the queue process instance
            other.updaters_queued.process_queue = Process_Queue(other=other)

    def init5(self,update=False):
        '''
        This model instance initialisation step is the last substantial one in which the dynamic solution of the DSGE
        model instance is finally computed using a choice of methods which can be called at runtime.
        
         
        :param other:     The DSGE model instance itother.
        :type other:      dsge_inst

        '''
        other = self._other
        txtpars = other.txtpars
        secs = txtpars.secs
        initlev = other._initlev
        ncpus = other._ncpus
        mk_hessian = other._mk_hessian
        mesg = other._mesg       
        #TODO: delay above and only import if needed
        from pymaclab.dsge.solvers.modsolvers import MODsolvers
        # Open the model solution tree branch
        other.modsolvers = MODsolvers()
        ################################### DYNARE++ ####################################
        # Open the dynare branch if dynare++ is installed on the system and in your PATH
        if dynare_flag:
            from pymaclab.dsge.solvers.modsolvers import Dynarepp
            other.modsolvers.dynarepp = Dynarepp({'_other':other,'_mesg':other._mesg})
        ############################# LINEAR METHODS !!! ################################
        # see if there are any log-linearized equations
        if all([False if 'None' in x else True for x in secs['modeq'][0]]):
            from pymaclab.dsge.solvers.modsolvers import PyUhlig, ForKlein
            if mesg: print "INIT: Computing DSGE model's log-linearized solution using Uhlig's Toolbox..."

            # Open the native Uhlig object
            intup = ((other.nendo,other.ncon,other.nexo,other.niid),
                 other.eqindx,
                 other.vreg,
                 other.llsys_list,
                 other.diffli1,
                 other.diffli2,
                 other._vtiming,
                 other.vardic)
            other.modsolvers.pyuhlig = PyUhlig(intup)

            # Open the Fortran Klein object
            intup = ((other.nendo,other.ncon,other.nexo,other.niid),
                 other.eqindx,
                 other.vreg,
                 other.llsys_list,
                 other.diffli1,
                 other.diffli2,
                 other._vtiming,
                 other.vardic)
            other.modsolvers.forklein = ForKlein(intup)
    ################## 1ST-ORDER NON-LINEAR METHODS !!! ##################
        if all([False if 'None' in x else True for x in secs['focs'][0]]):
            from pymaclab.dsge.solvers.modsolvers import ForKleinD
            
            if not update:
                if ncpus > 1 and mk_hessian:
                    if mesg: print "INIT: Computing DSGE model's Jacobian and Hessian using parallel approach..."
                    other.mkjahepp()
                elif ncpus > 1 and not mk_hessian:
                    if mesg: print "INIT: Computing DSGE model's Jacobian using parallel approach..."
                    other.mkjahepp()
                else:
                    if mesg: print "INIT: Computing DSGE model's Jacobian and Hessian using serial approach..."
                    other.mkjahe()
            else:
                if ncpus > 1 and mk_hessian:
                    if mesg: print "INIT: Computing DSGE model's Jacobian and Hessian using parallel approach..."
                    other.mkjaheppn()
                elif ncpus > 1 and not mk_hessian:
                    if mesg: print "INIT: Computing DSGE model's Jacobian using parallel approach..."
                    other.mkjaheppn()
                else:
                    if mesg: print "INIT: Computing DSGE model's Jacobian and Hessian using serial approach..."
                    other.mkjahen()                

            # Check if the obtained matrices A and B have correct dimensions
            if other.derivatives.jAA.shape[0] != other.derivatives.jAA.shape[1]:
                print "ERROR: Matrix A of derivatives does not have #vars=#equations"
            if other.derivatives.jBB.shape[0] != other.derivatives.jBB.shape[1]:
                print "ERROR: Matrix B of derivatives does not have #vars=#equations"

            # Open the Fortran KleinD object
            if 'nlsubsys' in dir(other):
                intup = (other.derivatives.numj,
                     other.nendo,other.nexo,
                     other.ncon,other.niid,other.sigma,
                     other.derivatives.jAA,other.derivatives.jBB,
                     other.vardic,other.vdic,
                     other.mod_name,other.audic,
                     other.derivatives.numjl,
                     other.nother)
            else:
                intup = (other.derivatives.numj,
                     other.nendo,other.nexo,
                     other.ncon,other.niid,other.sigma,
                     other.derivatives.jAA,other.derivatives.jBB,
                     other.vardic,other.vdic,
                     other.mod_name,other.audic)
            other.modsolvers.forkleind = ForKleinD(intup,other=other)
            # Make the AA and BB matrices as references available instead
            other.modsolvers.forkleind.A = other.derivatives.jAA
            other.modsolvers.forkleind.B = other.derivatives.jBB

    ################## 2ND-ORDER NON-LINEAR METHODS !!! ##################
        if all([False if 'None' in x else True for x in secs['vcvm'][0]]) and 'numh' in dir(other.derivatives):
            from pymaclab.dsge.solvers.modsolvers import PyKlein2D
            # Define intup
            if 'nlsubsys' in dir(other):
                intup = (other.derivatives.numj,other.derivatives.numh,
                     other.nendo,other.nexo,
                     other.ncon,other.niid,other.sigma,
                     other.derivatives.jAA,other.derivatives.jBB,
                     other.vardic,other.vdic,
                     other.mod_name,other.audic,
                     other.derivatives.numjl,other.derivatives.numhl,
                     other.nother)
            else:
                intup = (other.derivatives.numj,other.derivatives.numh,
                     other.nendo,other.nexo,
                     other.ncon,other.niid,other.sigma,
                     other.derivatives.jAA,other.derivatives.jBB,
                     other.vardic,other.vdic,
                     other.mod_name,other.audic)

            other.modsolvers.pyklein2d = PyKlein2D(intup,other=other)
            # Make the AA and BB matrices as references available instead
            other.modsolvers.pyklein2d.A = other.derivatives.jAA
            other.modsolvers.pyklein2d.B = other.derivatives.jBB
            other.modsolvers.pyklein2d.forkleind.A = other.derivatives.jAA
            other.modsolvers.pyklein2d.forkleind.B = other.derivatives.jBB
            
    def init_out(self):
        '''
        The final intializor section does some extra stuff after all has been done.
           
        :param self:     The DSGE model instance itself.
        :type self:      dsge_inst

        '''
        other = self._other
        initlev = other._initlev      
        # Make sure the jobserver has done his jobs
        jobserver.wait()