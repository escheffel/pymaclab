"""
GENERAL TEXTPARSER FOR MODFILES
"""
import os
import re
from ..tools import locate

modfpath = "" #TODO: get rid of these globals, make a config file

#TODO: make sure that everything gets evaluated in float terms
#for instance, in old cee.txt 1.03**(1/4) does int division
#TODO: rewrite to be a class factory and return a ParsedMod class, 
# because we don't need access to methods again.
class MODparser(object):
    """
    Parses a .mod file.

    Parameters
    ----------
    fname : str
        Filename to parse.
    modfmt : str {'pymaclab'}
        Only native .mod formats are currently supported.
    """
    def __init__(self, fname=None, modfmt='pymaclab'):
        self.type = type
        self.ffile = fname #TODO: why have both?
        self.filename = fname
        self.secs = {}     #TODO: change to sections
        input = open(os.path.join(modfpath, fname), 'r')
        wholefile = input.read()
        self.filestring = wholefile
        self.lines = self.filestring.splitlines()
        self.numlines = len(self.filestring.splitlines())

#Old section names
#        secnames = [['%Model Description','mod','MOD_loc'],
#                ['%Model Information','info','INF_loc'],
#                ['%Parameters','para','PA_loc'],
#                ['%Variable Vectors','varvec','VV_loc'],
#                ['%Boundary conditions','bocond','BC_loc'],
#                ['%Variable Substitution Non-Linear System','vsfocs',
#                    'VSFO_loc'],
#                ['%Non-linear first-order conditions','focs','FO_loc'],
#                ['%Steady States[closed-form]','sss','SS_loc'],
#                ['%Manual entry of sstate non-linear system','ssm','SSM_loc'],
#                ['%Log-Linearized Model Equations','modeq','ME_loc'],
#                ['%Variance-Covariance Matrix','vcvm','VCM_loc'],
#                ['%Minford Model Evaluation','mme','MME_loc']]

        secnames = [['%Model Information','info','INF_loc'],
                ['%Parameters','para','PA_loc'],
                ['%Variable Vectors','varvec','VV_loc'],
                ['%Boundary Conditions','bocond','BC_loc'],
                ['%Variable Substitution Non-Linear System','vsfocs',
                    'VSFO_loc'],
                ['%Non-Linear First-Order Conditions','focs','FO_loc'],
                ['%Steady States [Closed Form]','sss','SS_loc'],
                ['%Steady State Non-Linear System [Manual]','ssm','SSM_loc'],
                ['%Log-Linearized Model Equations','modeq','ME_loc'],
                ['%Variance-Covariance Matrix','vcvm','VCM_loc'],
                ['%Minford Model Evaluation','mme','MME_loc']]

        if modfmt == 'pymaclab':
            self.nativefmt(secnames)
        else:
            raise ValueError("modfmt %s not understood" % modfmt)

    def nativefmt(self, secnames):
        lines = self.lines
        numlines = self.numlines
        secs = self.secs        
        locdic = locate(lines,secnames)

        def readerfu(locname,locat,secname):
            list_tmp1 = []
            list_tmp2 = []
            row_iter=locat+1
            for x in lines[row_iter:]:
                x = re.sub("\s+", "", x)
                if '#' not in x and x != '' and x[0] != '%':
                    list_tmp1.append(lines[row_iter])
                    list_tmp2.append(lines[row_iter])
                    row_iter=row_iter+1
                elif '#' in x:
                    list_tmp2.append(lines[row_iter])
                    row_iter=row_iter+1
                elif x == '':
                    row_iter=row_iter+1
                elif x[0] == '%':
                    break
            while '' in list_tmp1:
                list_tmp1.remove('')
            self.secs[secname] = (list_tmp1,list_tmp2)

        for x in secnames:
            readerfu(x[1],locdic[x[1]],x[1])
