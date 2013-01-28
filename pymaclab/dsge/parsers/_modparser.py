'''
.. module:: _modparser
   :platform: Linux
   :synopsis: This is the (private) module responsible for collecting raw information from the DSGE model templated files. In here
              we don't make much use of Regex patterns as we are only extracting raw lines from the mod files.

.. moduleauthor:: Eric M. Scheffel <eric.scheffel@nottingham.edu.cn>


'''
import os
import re
from copy import deepcopy

modfpath = "" #TODO: get rid of these globals, make a config file

#TODO: make sure that everything gets evaluated in float terms
#for instance, in old cee.txt 1.03**(1/4) does int division
class ParsedMod(object):
    """
    Bunch pattern ParsedMod object

    Attributes
    ----------
    filename : str
        .mod file path and name
    secs : dict
        The sections of the .mod file and their contents
    filestring : str
        The .mod file as a string
    lines : list
        filestring split by lines
    numlines
        Length of line
    fmt: str
        The file format, i.e. pymaclab or dynare++
    """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

def locate(stringlines,varlist,regexm=False):
    '''
    A function, takes stringlines (split file) and a varlist
    and then creates a location dictionary for location of strings in file
    Can either use direct subset test or regex pattern matching test using regexm option=True
    '''
    locdic = {}
    for x in varlist:
        pattern = re.compile(x[0])
        row_iter=0
        while row_iter < len(stringlines):
            # not case sensitive
            if not regexm and x[0].lower() in stringlines[row_iter].lower():
                locdic[x[1]] = row_iter
                break
            elif regexm and pattern.search(stringlines[row_iter].lower()):
                locdic[x[1]] = row_iter
                break
            else:
                row_iter += 1
    return locdic


def read_nativefmt(filestring, fname):
    """
    Reads native .mod format and return as ParsedMod class.
    """
    secnames = [['%Model Information','info','INF_loc'],
                ['%Parameters','params','PA_loc'],
                ['%Variable Vectors','varvec','VV_loc'],
                ['%Boundary Conditions','bocond','BC_loc'],
                ['%Variable Substitution Non-Linear System','vsfocs','VSFO_loc'],
                ['%Non-Linear First-Order Conditions','focs','FO_loc'],
                ['%Steady States [Closed Form]','closedformss','SS_loc'],
                ['%Steady State Non-Linear System [Manual]','manualss','SSM_loc'],
                ['%Log-Linearized Model Equations','modeq','ME_loc'],
                ['%Variance-Covariance Matrix','vcvm','VCM_loc']]

    lines = filestring.splitlines()
    numlines = len(lines)
    secs = {}
    locdic = locate(lines,secnames)

    for x in secnames:
        section_name = x[1]
        list_tmp1 = []
        list_tmp2 = []
        row_iter = locdic[section_name]+1 # section starts 1 line after header
        for line in lines[row_iter:]:
            line = re.sub("\s+", "", line) # strip out any whitespace char
            if not line.startswith(('#','%')) and line != "":
                list_tmp1.append(lines[row_iter])
                list_tmp2.append(lines[row_iter])
                # really hold these duplicates if not labels?
            elif line.startswith('#'): # grab labels
                list_tmp2.append(lines[row_iter])
            elif line.startswith('%'):
                break
            row_iter += 1
        secs[section_name] = (list_tmp1,list_tmp2)
        
    #TODO: what do we actually _need_ and do we need a class?
    return ParsedMod(filename=fname, secs=secs, filestring=filestring, 
            lines=lines, numlines=len(lines),fmt='pymaclab')

def read_dynareppfmt(filestring, fname):
    """
    Reads dynarepp .mod format and return as ParsedMod class instance.
    """
    secnames = [['.*?var.*?;','var','VAR_loc'],
                ['.*?varexo.*?;','varexo','VAREXO_loc'],
                ['.*?parameters.*?;','parameters','PARAM_loc'],
                ['\w+.*?=.*?(?!\[)\d\.{0,1}\d*;','paramvals','PARAMVALS_loc'],
                ['.*?model.*?;','model','FOCS_loc'],
                ['.*?initval.*?;','initval','SSI_loc'],
                ['.*?vcov.*=?','vcov','VCOV_loc'],
                ['.*?order.*?=.*?\d{1,1}','order','ORDER_loc']]

    lines = filestring.splitlines()
    numlines = len(lines)
    secs = {}
    locdic = locate(lines,secnames,regexm=True)

    for i1,x in enumerate(secnames):
        section_name = x[1]
        list_tmp1 = []
        list_tmp2 = []
        # section starts 1 line after header, but only for "model" and "initval" sections!
        if x[1] == 'model' or x[1] == 'initval': row_iter = deepcopy(locdic[section_name]+1)
        else: row_iter = deepcopy(locdic[section_name])
        for line in lines[row_iter:]:
            if x[1] == 'order' and line != '':
                list_tmp1.append(lines[row_iter])
                list_tmp2.append(lines[row_iter])
            elif x[1] == 'parameters' and ';' not in line:
                list_tmp1.append(lines[row_iter])
                list_tmp2.append(lines[row_iter])
            elif x[1] == 'parameters' and ';' in line:
                list_tmp1.append(lines[row_iter])
                list_tmp2.append(lines[row_iter])                
                break
            elif x[1] == 'varexo' and ';' not in line:
                list_tmp1.append(lines[row_iter])
                list_tmp2.append(lines[row_iter])
            elif x[1] == 'varexo' and ';' in line:
                list_tmp1.append(lines[row_iter])
                list_tmp2.append(lines[row_iter])                
                break
            elif x[1] == 'var' and ';' not in line:
                list_tmp1.append(lines[row_iter])
                list_tmp2.append(lines[row_iter])
            elif x[1] == 'var' and ';' in line:
                list_tmp1.append(lines[row_iter])
                list_tmp2.append(lines[row_iter])                
                break
            elif x[1] == 'paramvals' and line != '':
                list_tmp1.append(lines[row_iter])
                list_tmp2.append(lines[row_iter])
            elif x[1] == 'paramvals' and line == '':
                break
            elif x[1] == 'model' and 'end;' not in line:
                list_tmp1.append(lines[row_iter])
                list_tmp2.append(lines[row_iter])
            elif x[1] == 'model' and 'end;' in line:
                break
            elif x[1] == 'initval' and 'end;' not in line:
                list_tmp1.append(lines[row_iter])
                list_tmp2.append(lines[row_iter])
            elif x[1] == 'initval' and 'end;' in line:
                break
            elif x[1] == 'vcov' and '];' not in line:
                list_tmp1.append(lines[row_iter])
                list_tmp2.append(lines[row_iter])
            elif x[1] == 'vcov' and '];' in line:
                list_tmp1.append(lines[row_iter])
                list_tmp2.append(lines[row_iter])                
                break
            row_iter += 1
        secs[section_name] = (list_tmp1,list_tmp2)
        
    #TODO: what do we actually _need_ and do we need a class?
    return ParsedMod(filename=fname, secs=secs, filestring=filestring, 
            lines=lines, numlines=len(lines),fmt='dynarepp')

            

def parse_mod(fname, modfmt='pymaclab'):
#class MODparser(object):
    """
    Parses a .mod file.

    Parameters
    ----------
    fname : str
        Filename to parse. Could also be the model info as string itself.
    modfmt : str {'pymaclab'}
        Only native .mod formats are currently supported.

    Returns
    -------
    ParsedMod class.
    """
    # Need to check here what kind of file has been read in, pymaclab or dynare++ or something else?
    # This is a pattern which should identify dynare++ model file
    dynpat = re.compile('.*?end;') 
    if os.path.exists(fname):
        with open(os.path.join(modfpath, fname), 'r') as txtfile:
            filestring = txtfile.read()           
            if dynpat.search(filestring): modfmt = 'dynarepp'
            if modfmt == 'dynarepp':
                parsed_mod = read_dynareppfmt(filestring, fname)
            elif modfmt == 'pymaclab':
                parsed_mod = read_nativefmt(filestring, fname)
            else:
                raise ValueError("modfmt %s not understood" % modfmt)
    # Assume here that the model information was directly passed as a string
    elif '\n' in fname and len(fname) > 50:
        if dynpat.search(fname): modfmt = 'dynarepp'
        if modfmt == 'pymaclab':
            parsed_mod = read_nativefmt(fname, fname)
        elif modfmt == 'dynarepp':
            parsed_mod = read_dynareppfmt(fname, fname)
        else:
            raise ValueError("modfmt %s not understood" % modfmt)
    #  Otherwise don't know what to do with the input and throw error
    else:
        raise ValueError("modfmt %s not understood" % modfmt)
        
    return parsed_mod
