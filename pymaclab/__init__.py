'''
.. module:: pymaclab
   :platform: Linux
   :synopsis: The var module contains the VAR class for estimating and doing further work with Vector Autoregressions commonly
              used in applied macroeconometrics. It supports advanced methods such as bootstrapping confidence intervals including
              Killian's boostrap-after-bootstrap small-sample bias correction. Also CPU-intensive methods such as the bootstrap can
              be computed using Parallel Python to exploit multi-core CPUs. Pretty plotting methods are also included which depend
              on matplotlib.

.. moduleauthor:: Eric M. Scheffel <eric.scheffel@nottingham.edu.cn>


'''

import os as OPS
from dsge import macrolab
from stats import var
from stats import favar
import linalg
import sys

# Expose the version number into library root
from version import version as __version__

# Some helper functions
def db_graph(dbase,tseries):
	fig = TPL.tsfigure()
	fsp = fig.add_tsplot(111)
	P.title(dbase.datdic[tseries]['Dat_Desc'])
	fsp.tsplot(dbase.datdic[tseries]['infile'], '-', label=tseries)
	fsp.format_dateaxis()
	fsp.set_xlim(int(dbase.datdic[tseries]['infile'].start_date),
		     int(dbase.datdic[tseries]['infile'].end_date))
	fsp.legend()
	P.show()

def newMOD(txtfile=None,dbase=None,initlev=2,mesg=False,ncpus=1,mk_hessian=True,use_focs=False,ssidic=None):
	'''
	Model's second intialisation method called by newMOD() function call. The model's
	__init__() method only creates the instance and adds information, but does no
	parsing and computations whatsoever. init2() parses and ready's the model for
	calculating the steady state and the dynamic solution to the model.
	
	init_lev = 0: only parsing, no calculations, but preparing for manual SS solution
	init_lev = 1: parsing and steady state calculations
	init_lev = 2: parsing, steady state calculations and dynamic solution computation
	'''
	modobj = macrolab.DSGEmodel(txtfile,dbase=dbase,initlev=initlev,mesg=mesg,ncpus=ncpus,
	                            mk_hessian=mk_hessian,use_focs=use_focs,ssidic=ssidic)
	modobj.init1()
	modobj.init1a()
	modobj.init1b()
	modobj.init1c()
	modobj.init2()
	# Ready to solve for SS manually at command prompt
	if initlev == 0:
		modobj.init_out()
		return modobj
	# SS solved automatically
	modobj.init3()
	if initlev == 1:
		modobj.init_out()
		return modobj
	# Dynamic system prepared and solved
	modobj.init4()
	modobj.init5()
	if initlev == 2:
		modobj.init_out()
		return modobj

def newDB():
	return macrolab.TSDataBase()

def newVAR(data=None,vnames=None,pnames=None,svnames=None,irfs=True,boot=True,plot=True,conf=None,mesg=False):
	return var.VAR(data=data,vnames=vnames,pnames=pnames,svnames=svnames,irfs=irfs,boot=boot,plot=plot,conf=conf,mesg=mesg)

def newFAVAR(dates=None,rdates=None,data=None,freq='M',vnames=None,pnames=None,svnames=None,irfs=True,rescale=False,boot=True,plot=True,sfacs='auto',init=None,conf=None,mesg=False):
	return favar.FAVAR(dates=dates,rdates=rdates,data=data,freq=freq,vnames=vnames,pnames=pnames,svnames=svnames,irfs=irfs,rescale=rescale,boot=boot,plot=plot,sfacs=sfacs,init=init,conf=conf,mesg=mesg)

def modinfo(model):
	offset = 40
	# Print the Description
	print '\n'
	print'Model Description:'
	print 60*'-'
	for x in model.txtpars.secs['mod'][1]:
		print x
	print 60*'-'
	# What's it's name?
	print 'Model Name:'+(offset-len('Model Name:'))*' '+model.modname
	# From which File parsed?
	print 'Created from:'+(offset-len('Created from:'))*' '+model.modfile
	# Supplied with data?
	if 'data' in dir(model) and model.data != None:
		print 'Data attached?'+(offset-len('Data attached?'))*' '+'YES'
	else:
		print 'Data attached?'+(offset-len('Data attached?'))*' '+'NO'
	# Possesses matuhlig?
	if 'matuhlig' in dir(model.modsolvers):
		# Possesses matuhlig solution?
		if 'PP' in dir(model.modsolvers.matuhlig):
			print 'MatUhlig Solution available?'+(offset-len('MatUhlig Solution available?'))*' '+'YES(*)'
		else:
			print 'MatUhlig Solution available?'+(offset-len('MatUhlig Solution available?'))*' '+'NO(*)'
	else:
		print 'MatUhlig Solution available?'+(offset-len('MatUhlig Solution available?'))*' '+'NO(X)'
	# Possesses pyuhlig?
	if 'pyuhlig' in dir(model.modsolvers):
		# Possesses pyuhlig solution?
		if 'PP' in dir(model.modsolvers.pyuhlig):
			print 'PyUhlig Solution available?'+(offset-len('PyUhlig Solution available?'))*' '+'YES(*)'
		else:
			print 'PyUhlig Solution available?'+(offset-len('PyUhlig Solution available?'))*' '+'NO(*)'
	else:
		print 'PyUhlig Solution available?'+(offset-len('PyUhlig Solution available?'))*' '+'NO(X)'
	# Possesses matklein?
	if 'matklein' in dir(model.modsolvers):
		# Possesses matklein solution?
		if 'P' in dir(model.modsolvers.matklein):
			print 'MatKlein Solution available?'+(offset-len('MatKlein Solution available?'))*' '+'YES(*)'
		else:
			print 'MatKlein Solution available?'+(offset-len('MatKlein Solution available?'))*' '+'NO(*)'
	else:
		print 'MatKlein Solution available?'+(offset-len('MatKlein Solution available?'))*' '+'NO(X)'
	# Possesses forklein?
	if 'forklein' in dir(model.modsolvers):
		# Possesses forklein solution?
		if 'P' in dir(model.modsolvers.forklein):
			print 'ForKlein Solution available?'+(offset-len('ForKlein Solution available?'))*' '+'YES(*)'
		else:
			print 'ForKlein Solution available?'+(offset-len('ForKlein Solution available?'))*' '+'NO(*)'
	else:
		print 'ForKlein Solution available?'+(offset-len('ForKlein Solution available?'))*' '+'NO(X)'
	# Possesses forkleind?
	if 'forkleind' in dir(model.modsolvers):
		# Possesses forkleind solution?
		if 'P' in dir(model.modsolvers.forkleind):
			print 'ForKleinD Solution available?'+(offset-len('ForKleinD Solution available?'))*' '+'YES(*)'
		else:
			print 'ForKleinD Solution available?'+(offset-len('ForKleinD Solution available?'))*' '+'NO(*)'
	else:
		print 'ForKleinD Solution available?'+(offset-len('ForKleinD Solution available?'))*' '+'NO(X)'
	# Possesses matklein2d?
	if 'matklein2d' in dir(model.modsolvers):
		# Possesses matklein2d solution?
		if 'P' in dir(model.modsolvers.matklein2d):
			print 'MatKlein2D Solution available?'+(offset-len('MatKlein2D Solution available?'))*' '+'YES(*)'
		else:
			print 'MatKlein2D Solution available?'+(offset-len('MatKlein2D Solution available?'))*' '+'NO(*)'
	else:
		print 'MatKlein2D Solution available?'+(offset-len('MatKlein2D Solution available?'))*' '+'NO(X)'

def modsolve(model,stype):
	if stype == 'matuhlig':
		if 'PP' in dir(model.modsolvers.matuhlig):
			del model.modsolvers.matuhlig.PP
			del model.modsolvers.matuhlig.QQ
			del model.modsolvers.matuhlig.RR
			del model.modsolvers.matuhlig.SS
			del model.modsolvers.matuhlig.WW
		model.modsolvers.matuhlig.solve()
		if 'PP' in dir(model.modsolvers.matuhlig):
			return dict(PP=model.modsolvers.matuhlig.PP,
				    QQ=model.modsolvers.matuhlig.QQ,
				    RR=model.modsolvers.matuhlig.RR,
				    SS=model.modsolvers.matuhlig.SS,
				    WW=model.modsolvers.matuhlig.WW)
		else:
			print 'Error: Solution not obtained!'
			return False
	elif stype == 'matklein':
		if 'P' in dir(model.modsolvers.matklein):
			del model.modsolvers.matklein.P
			del model.modsolvers.matklein.F
		model.modsolvers.matklein.solve()
		if 'P' in dir(model.modsolvers.matklein):
			return dict(P=model.modsolvers.matklein.P,
				    F=model.modsolvers.matklein.F)
		else:
			print 'Error: Solution not obtained!'
			return False
	elif stype == 'pyuhlig':
		if 'PP' in dir(model.modsolvers.pyuhlig):
			del model.modsolvers.pyuhlig.PP
			del model.modsolvers.pyuhlig.QQ
			del model.modsolvers.pyuhlig.RR
			del model.modsolvers.pyuhlig.SS
			del model.modsolvers.pyuhlig.WW
		model.modsolvers.pyuhlig.solve()
		if 'PP' in dir(model.modsolvers.pyuhlig):
			return dict(PP=model.modsolvers.pyuhlig.PP,
				    QQ=model.modsolvers.pyuhlig.QQ,
				    RR=model.modsolvers.pyuhlig.RR,
				    SS=model.modsolvers.pyuhlig.SS,
				    WW=model.modsolvers.pyuhlig.WW)
		else:
			print 'Error: Solution not obtained!'
			return False
	elif stype == 'forklein':
		if 'P' in dir(model.modsolvers.forklein):
			del model.modsolvers.forklein.P
			del model.modsolvers.forklein.F
		model.modsolvers.forklein.solve()
		if 'P' in dir(model.modsolvers.forklein):
			return dict(P=model.modsolvers.forklein.P,
				    F=model.modsolvers.forklein.F)
		else:
			print 'Error: Solution not obtained!'
			return False
	else:
		return 'Solution method not recognized!'

def lmods():
	aa = [x[0] for x in filter(lambda x: isinstance(x[1],macrolab.DSGEmodel),macrolab.locdic.items())]
	if len(aa) == 0:
		return 'You have no maclab.DSGE models in memory yet!'
	aa.sort()
	tabu = 5
	# Find longest model name
	mlen = 0
	for x in aa:
		if len(x) > mlen:
			mlen = len(x)
	# Find longest model description
	mdlen = 0
	for x in aa:
		if len(macrolab.locdic[x].modname) > mdlen:
			mdlen = len(macrolab.locdic[x].modname)
	# Find longest model var number
	mvlen = 0
	for x in aa:
		if len(str(macrolab.locdic[x].nall)) > mvlen:
			mvlen = len(str(macrolab.locdic[x].nall))
	# Which is bigger? #Vars or mvlen?
	if len('#Vars') > mvlen:
		mvlen = len('#Vars')

	print ''
	print 'Name'+(mlen+tabu-len('Name'))*' '+'Description'+(mdlen+tabu-len('Description'))*' '+'#Vars'
	print '-'*(mlen+mdlen+mvlen+2*tabu)
	for x in aa:
		print x+(mlen+tabu-len(x))*' '+macrolab.locdic[x].modname+(mdlen+tabu-len(macrolab.locdic[x].modname))*' '+str(macrolab.locdic[x].nall)

def ldbs():
	aa = [x[0] for x in filter(lambda x: isinstance(x[1],macrolab.TSDataBase),macrolab.locdic.items())]
	if len(aa) == 0:
		return 'You have no maclab.DB databases in memory yet!'
	aa.sort()
	return aa

def lvars():
	aa = [x[0] for x in filter(lambda x: isinstance(x[1],macrolab.VAR),macrolab.locdic.items())]
	if len(aa) == 0:
		return 'You have no maclab.VAR vector autoregressions in memory yet!'
	aa.sort()
	return aa

def pyed(pyfile=None):
	cwd = OPS.getcwd()
	dirli = OPS.listdir(cwd)
	if pyfile == None:
		return 'You did not specify the Python file to edit!'
	elif pyfile+'.py' in dirli:
		cmd = macrolab.pyedpath+' '+pyfile+'.py'
		OPS.system(cmd)
	else:
		return 'The file '+pyfile+'.py '+'was not found in current directory!'

def modedit(model):
	model.txted()

def texedit(model):
	model.texed()

def explain(model):
	model.pdf()


# Import this late because cascading imports need newMOD() function
import modfiles
from .modfiles.makemod import make_modfile
