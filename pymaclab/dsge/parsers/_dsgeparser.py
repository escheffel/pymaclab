'''
.. module:: _dsgeparser
   :platform: Linux
   :synopsis: This is the (private) module responsible for carefully extracting meaningful information from the DSGE model templated
              files. In here we make much use of Regex patterns as we are extracting information from the mod file lines.

.. moduleauthor:: Eric M. Scheffel <eric.scheffel@nottingham.edu.cn>


'''
import re
from copy import deepcopy
import numpy.matlib as mat
import sys
import os
# This used to import sympycore, but should now also work with sympy
import sympycore as SP


#TODO: why have these as nested functions?
def populate_model_stage_one(self, secs):
    """
    1st stage population of DSGE model.  Does not need Steady State.
    """
    
    # This is a special dictionary which can be handed over to the template engines (i.e. Jinja2)
    if 'template_paramdic' not in dir(self): self.template_paramdic = {}


    # Get all information from the varvec section
    if all([False if 'None' in x else True for x in secs['varvec'][0]]):
        # match something like [1]  z(t):eps_z(t):techshock{exo}[log,hp]
        vtexp =  re.compile('^.*?\[.*?\]\s*(?P<vari>.+?)\s*:\s*(?P<varn>.+?)\s*\{(?P<vtype>[^\[]+)}\s*(?P<mod>\[.*?\]){0,1}')
        # match something like [15] x2(t):wrec2{con}[log,hp] OR [1] z(t):eps_z(t):techshock{exo}[log,hp]
        viiexp = re.compile('^.*?\[.*?\]\s*(?P<vari>.+?)\s*:\s*(?P<viid>.+?)\s*:\s*(?P<varn>.+?)\s*\{(?P<vtype>[^\[]+)}\s*(?P<mod>\[.*?\]){0,1}')
        voexp =  re.compile('^.*?\[.*?\]\s*(?P<vari>@.+?):(?P<varn>[^\[]+)\s*(?P<mod>\[.*?\]){0,1}')
        vardic = {}
        audic = {}
        vardic['endo'] = {}
        audic['endo'] = {}
        vardic['exo'] = {}
        audic['exo'] = {}
        vardic['con'] = {}
        audic['con'] = {}
        vardic['other'] = {}
        vardic['endo']['var'] = []
        vardic['endo']['mod'] = []
        vardic['exo']['var'] = []
        vardic['exo']['mod'] = []
        vardic['con']['var'] = []
        vardic['con']['mod'] = []
        vardic['other']['var'] = []
        vardic['other']['mod'] = []
        audic['endo']['var'] = []
        audic['endo']['mod'] = []
        audic['con']['var'] = []
        audic['con']['mod'] = []
        audic['exo']['var'] = []
        audic['exo']['mod'] = []

        for x in secs['varvec'][0]:
            if viiexp.search(x):
                ma = viiexp.search(x)
                vari = ma.group('vari').strip()
                viid = ma.group('viid').strip()
                varn = ma.group('varn').strip()
                vtype = ma.group('vtype').strip()
                mods = ma.group('mod')
                vardic[vtype]['var'].append([vari,varn,viid])
                if mods != None:
                    mods = mods.strip()
                    if ',' in mods:
                        vardic[vtype]['mod'].append(mods[1:-1].strip().split(','))
                    else:
                        vardic[vtype]['mod'].append([mods[1:-1].strip(),])
                else:
                    vardic[vtype]['mod'].append([])

            elif vtexp.search(x):
                ma = vtexp.search(x)
                vari = ma.group('vari').strip()
                varn = ma.group('varn').strip()
                vtype = ma.group('vtype').strip()
                mods = ma.group('mod')
                vardic[vtype]['var'].append([vari,varn])
                if mods != None:
                    mods = mods.strip()
                    if ',' in mods:
                        vardic[vtype]['mod'].append(mods[1:-1].strip().split(','))
                    else:
                        vardic[vtype]['mod'].append([mods[1:-1].strip(),])
                else:
                    vardic[vtype]['mod'].append([])

            elif voexp.search(x):
                ma = voexp.search(x)
                # Slice off the @
                vari = ma.group('vari')[1:].strip()
                varn = ma.group('varn').strip()
                mods = ma.group('mod')
                vardic['other']['var'].append([vari,varn])
                if mods != None:
                    mods = mods.strip()
                    if ',' in mods:
                        vardic['other']['mod'].append(mods[1:-1].strip().split(','))
                    else:
                        vardic['other']['mod'].append([mods[1:-1].strip(),])
                else:
                    vardic['other']['mod'].append([])

        self.nendo = len(vardic['endo']['var'])
        self.nexo = len(vardic['exo']['var'])
        self.ncon = len(vardic['con']['var'])
        self.nother = len(vardic['other']['var'])
        self.nstat = self.nendo+self.nexo
        self.nall = self.nstat+self.ncon
        self.vardic = vardic
        self.audic = audic
        # Save for template instantiation
        self.template_paramdic['vardic'] = deepcopy(vardic)
    else:
        self.template_paramdic['vardic'] = False

    # Extract the model name and description
    if all([False if 'None' in x else True for x in secs['info'][0]]):
        for x in secs['info'][0]:
            if 'Name' in x:
                self.mod_name = x.split('=')[1].replace(';','').strip()
                # Save for template instantiation
                self.template_paramdic['mod_name'] = deepcopy(self.mod_name)          
            if 'Desc' in x:
                self.mod_desc = x.split('=')[1].replace(';','').strip()
                # Save for template instantiation
                self.template_paramdic['mod_desc'] = deepcopy(self.mod_desc)
        if not self.template_paramdic.has_key('mod_name'):
            self.template_paramdic['mod_name'] = False
        if not self.template_paramdic.has_key('mod_desc'):
            self.template_paramdic['mod_desc'] = False

    # Extract parameters into dictionary
    if all([False if 'None' in x else True for x in secs['params'][0]]):
        param = {}
        # need to do this so users don't have to worry about integer division
        # but preserves integer division for sympycore stuff
        safe_div = """
from __future__ import division
"""
        for x in secs['params'][0]:
            if ']' in x: x = x.split(']')[1].lstrip()
            list_tmp = x.split(';')
            list_tmp = list_tmp[0].split('=')[:]
            str_tmp1 = list_tmp[0].strip()
            str_tmp2 = list_tmp[1].strip()
            #NOTE: this *should* be safe, but users should know what's
            # in the .mod file
            exec(safe_div + "param['"+str_tmp1+"']=" + str_tmp2, {}, locals())
            locals()[str_tmp1] = param[str_tmp1]
        self.paramdic = param
        # Save for template instantiation
        self.template_paramdic['paramdic'] = deepcopy(param)
    else:
        self.template_paramdic['paramdic'] = False

    # Collect information on manual (closed-form) steady state
    if all([False if 'None' in x else True for x in secs['closedformss'][0]]):
        # Join multiline steady state definitions
        mansys = secs['closedformss'][0]
        list_tmp1 = []
        i1=0
        counter=0
        for x in mansys:
            if '...' in x:
                counter = counter + 1
            elif '...' not in x:
                if counter == 0:
                    # Just in case people have used [1] markers...
                    if ']' in x: x = x.split(']')[1].lstrip().rstrip()
                    list_tmp1.append(x)
                elif counter > 0:
                    str_tmp = ''
                    for y in mansys[i1-counter:i1+1]:
                        str_tmp = str_tmp + y.replace('...','')
                    # Just in case people have used [1] markers...
                    if ']' in str_tmp: str_tmp = str_tmp.split(']')[1].lstrip().rstrip()
                    list_tmp1.append(str_tmp)
                    counter = 0 
            i1=i1+1
        self.manss_sys = list_tmp1
        # Save for template instantiation
        self.template_paramdic['manss_sys'] = deepcopy(list_tmp1)
    else:
        self.template_paramdic['manss_sys'] = False
        
    # Collect info on numerical steady state
    if all([False if 'None' in x else True for x in secs['manualss'][0]]):
        _mreg_alt = '^\s*[a-zA-Z]*_bar(?!=\+|-|\*|/])\s*=\s*.*'
        _mreg = '\[\d+\]\s*[a-zA-Z]*_bar(?!=\+|-|\*|/])\s*=\s*.*'
        _mreg2 = '\[\d+\]\s*[a-zA-Z]*(?!=\+|-|\*|/])\s*=\s*[0-9]*\.[0-9]*'
        _mregfocs = 'USE_FOCS=\[.*?\]'
        mregfocs = re.compile(_mregfocs)
        use_focs = False
        for lino in secs['manualss'][0]:
            if mregfocs.search(lino):
                # Save as the list of equations the be used
                use_focs = eval(lino.split('=')[1].replace(';','').strip())
                self._internal_focs_used = True
        mreg = re.compile(_mreg_alt+'|'+_mreg+'|'+_mreg2)
        indx = []
        ssidic={}
        list_tmp = []
        anyssi_found = False
        i1=0
        counter=0
        for x in secs['manualss'][0]:
            if mreg.search(x):
                raw_str = x.replace(';','')
                anyssi_found = True
                ma = mreg.search(raw_str)
                if ']' in ma.group(): str1 = ma.group().replace(';','').split('=')[0].split(']')[1].strip()
                else: str1 = ma.group().replace(';','').split('=')[0].lstrip().rstrip()
                str2 = ma.group().split('=')[1].strip()
                ssidic[str1] = eval(str2)
                # Expose the evaluated values for recursive smart evaluation
                locals().update(self.paramdic)
                locals()[str1] = eval(str2)
                indx.append(i1)
            elif not mreg.search(x) and '...' in x:
                counter = counter + 1
            elif not mreg.search(x) and '...' not in x:
                if counter == 0:
                    if ']' in x: list_tmp.append(x.replace(';','').split(']')[1].split('=')[0].strip())
                    else: list_tmp.append(x.replace(';','').split('=')[0].strip())
                elif counter > 0:
                    str_tmp = ''
                    for y in secs['manualss'][0][i1-counter:i1+1]:
                        str_tmp = str_tmp + y.replace('...','').strip()
                    if ']' in str_tmp: list_tmp.append(str_tmp.split(']')[1].split('=')[0].replace(';','').strip())
                    else: list_tmp.append(str_tmp.split('=')[0].replace(';','').strip())
                    counter = 0 
            i1=i1+1
        if anyssi_found:
            self.ssidic = deepcopy(ssidic)
            # Save for template instantiation
            self.template_paramdic['ssidic'] = deepcopy(ssidic)            
        else:
            #Make empty to be filled by other procedure
            self.ssidic = {}
            # Save for template instantiation
            self.template_paramdic['ssidic'] = False

        if not use_focs:
            self.ssys_list = deepcopy(list_tmp) 
            # Save for template instantiation
            self.template_paramdic['ssys_list'] = deepcopy(list_tmp)
            self.template_paramdic['use_focs'] = False
        elif use_focs and anyssi_found:
            self._use_focs = deepcopy(use_focs)
            self._ssidic = deepcopy(ssidic)
            # Save for template instantiation
            self.template_paramdic['ssys_list'] = False
            self.template_paramdic['use_focs'] = deepcopy(use_focs)
        elif use_focs and not anyssi_found:
            self._use_focs = deepcopy(use_focs)
            # Save for template instantiation
            self.template_paramdic['ssys_list'] = False
            self.template_paramdic['use_focs'] = use_focs
            
    else:
        # Save for template instantiation
        self.template_paramdic['ssidic'] = False       
        self.template_paramdic['ssys_list'] = False
        self.template_paramdic['use_focs'] = False
        

    return self


############# BELOW HERE IS ALL FOR 2ND STAGE ###########

def mkaug1(self, insys,othersys):         
    # Determine the lengths of augmented vars

    list_tmp2 = deepcopy(insys)
    list_tmp1 = deepcopy(othersys)

    endosp = []
    for x in self.vardic['endo']['var']:
        endosp = endosp + [[x,[-1,0]]]
    exosp = []
    for x in self.vardic['exo']['var']:
        exosp = exosp + [[x,[0,1]]]
    consp = []
    for x in self.vardic['con']['var']:
        consp = consp + [[x,[0,1]]]

    alldic = {}
    alldic.update(self.paramdic)
    alldic.update(self.sstate)

    spvdic = {}
    spvdic['endo'] = endosp
    spvdic['exo'] = exosp
    spvdic['con'] = consp

    spvdic2 = deepcopy(spvdic)

    endoli = [x[0].split('(')[0].strip() for x in self.vardic['endo']['var']]
    exoli = [x[0].split('(')[0].strip() for x in self.vardic['exo']['var']]
    conli = [x[0].split('(')[0].strip() for x in self.vardic['con']['var']]


    patup = ('{-100,100}|None','all','{-100,100}')
    count = 0
    for i1,line in enumerate(list_tmp2):
        iterob = self.vreg(patup,line,True,'max')
        if iterob:
            iterob = list(iterob)
        else:
            continue
        iterob.reverse()
        for vma in iterob:
            vtype = vma[1][0]
            vartime = vma[2][2]
            vari = vma[2][1]
            pos = vma[3][0]
            poe = vma[3][1]
            if vtype == 'endo':
                indx = endoli.index(vari)
            elif vtype == 'con':
                indx = conli.index(vari)
            elif vtype == 'exo':
                indx = exoli.index(vari)
            else:
                continue

            # Check for endo
            if vtype == 'endo' and int(vartime) > spvdic['endo'][indx][1][1]:
                spvdic2['endo'][indx][1][1] = int(vartime)
                tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
                newvar = vari+'_F'+tind+'(t)'
                newname =  vari+'_F'+str(abs(int(vartime)))
                list_tmp2[i1] = list_tmp2[i1][:pos]+newvar.split('(')[0]+'(t)'+list_tmp2[i1][poe:]
                for i2 in range(int(vartime)):
                    tind = (5-len(str(abs(i2+1)-1)))*'0'+str(abs(i2+1)-1)
                    newvar = vari+'_F'+tind+'(t)'
                    newname =  vari+'_F'+str(abs(i2))
                    if [newvar,newname] not in self.vardic['con']['var']:
                        self.vardic['con']['var'].append([newvar,newname])
                        self.vardic['con']['mod'].append(self.vardic['endo']['mod'][indx])
                        self.audic['con']['var'].append([newvar,newname])
                        self.audic['con']['mod'].append(self.vardic['endo']['mod'][indx])
                        if self.sstate.has_key(vari+'_bar'):
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
                        else:
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
                    continue
            elif vtype == 'endo' and int(vartime) < spvdic['endo'][indx][1][0]:
                spvdic2['endo'][indx][1][0] = int(vartime)
                tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
                newvar = vari+'_B'+tind+'(t)'
                newname =  vari+'_B'+str(abs(int(vartime)))
                list_tmp2[i1] = list_tmp2[i1][:pos]+newvar.split('(')[0]+'(t-1)'+list_tmp2[i1][poe:]
                for i2 in range(1,abs(int(vartime))):
                    tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
                    newvar = vari+'_B'+tind+'(t)'
                    newname =  vari+'_B'+str(abs(i2+1))
                    if [newvar,newname] not in self.vardic['endo']['var']:
                        self.vardic['endo']['var'].append([newvar,newname])
                        self.vardic['endo']['mod'].append(self.vardic['endo']['mod'][indx])
                        self.audic['endo']['var'].append([newvar,newname])
                        self.audic['endo']['mod'].append(self.vardic['endo']['mod'][indx])
                        if self.sstate.has_key(vari+'_bar'):
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
                        else:
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
                    continue
            # Check for exo
            if vtype == 'exo' and int(vartime) > spvdic['exo'][indx][1][1]:
                spvdic2['exo'][indx][1][1] = int(vartime)
                tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
                newvar = vari+'_F'+tind+'(t)'
                newname =  vari+'_F'+str(abs(int(vartime)-1))
                list_tmp2[i1] = list_tmp2[i1][:pos]+newvar.split('(')[0]+'(t)'+list_tmp2[i1][poe:]
                for i2 in range(int(vartime)):
                    tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
                    newvar = vari+'_F'+tind+'(t)'
                    newname =  vari+'_F'+str(abs(i2+1))
                    if [newvar,newname] not in self.vardic['con']['var']:
                        self.vardic['con']['var'].append([newvar,newname])
                        self.vardic['con']['mod'].append(self.vardic['exo']['mod'][indx])
                        self.audic['con']['var'].append([newvar,newname])
                        self.audic['con']['mod'].append(self.vardic['exo']['mod'][indx])
                        if self.sstate.has_key(vari+'_bar'):
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
                        else:
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
                    continue
            elif vtype == 'exo' and int(vartime) < spvdic['exo'][indx][1][0]:
                spvdic2['exo'][indx][1][0] = int(vartime)
                tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
                newvar = vari+'_B'+tind+'(t)'
                newname =  vari+'_B'+str(abs(int(vartime)))
                list_tmp2[i1] = list_tmp2[i1][:pos]+newvar.split('(')[0]+'(t-1)'+list_tmp2[i1][poe:]
                for i2 in range(abs(int(vartime))):
                    tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
                    newvar = vari+'_B'+tind+'(t)'
                    newname =  vari+'_B'+str(abs(i2+1))
                    if [newvar,newname] not in self.vardic['endo']['var']:
                        self.vardic['endo']['var'].append([newvar,newname])
                        self.vardic['endo']['mod'].append(self.vardic['exo']['mod'][indx])
                        self.audic['endo']['var'].append([newvar,newname])
                        self.audic['endo']['mod'].append(self.vardic['exo']['mod'][indx])
                        if self.sstate.has_key(vari+'_bar'):
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
                        else:
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
                    continue
            # Check for con
            if vtype == 'con' and int(vartime) > spvdic['con'][indx][1][1]:
                spvdic2['con'][indx][1][1] = int(vartime)
                tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
                newvar = vari+'_F'+tind+'(t)'
                newname =  vari+'_F'+str(abs(int(vartime)-1))
                list_tmp2[i1] = list_tmp2[i1][:pos]+newvar.split('(')[0]+'(t)'+list_tmp2[i1][poe:]
                for i2 in range(int(vartime)):
                    tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
                    newvar = vari+'_F'+tind+'(t)'
                    newname =  vari+'_F'+str(abs(i2+1))
                    if [newvar,newname] not in self.vardic['con']['var']:
                        self.vardic['con']['var'].append([newvar,newname])
                        self.vardic['con']['mod'].append(self.vardic['con']['mod'][indx])
                        self.audic['con']['var'].append([newvar,newname])
                        self.audic['con']['mod'].append(self.vardic['con']['mod'][indx])
                        if self.sstate.has_key(vari+'_bar'):
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
                        else:
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
                    continue
            elif vtype == 'con' and int(vartime) < spvdic['con'][indx][1][0]:
                spvdic2['con'][indx][1][0] = int(vartime)
                tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
                newvar = vari+'_B'+tind+'(t)'
                newname =  vari+'_B'+str(abs(int(vartime)))
                list_tmp2[i1] = list_tmp2[i1][:pos]+newvar.split('(')[0]+'(t-1)'+list_tmp2[i1][poe:]
                for i2 in range(abs(int(vartime))):
                    tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
                    newvar = vari+'_B'+tind+'(t)'
                    newname =  vari+'_B'+str(abs(i2+1))
                    if [newvar,newname] not in self.vardic['endo']['var']:
                        self.vardic['endo']['var'].append([newvar,newname])
                        self.vardic['endo']['mod'].append(self.vardic['con']['mod'][indx])
                        self.audic['endo']['var'].append([newvar,newname])
                        self.audic['endo']['mod'].append(self.vardic['con']['mod'][indx])
                        self.sstate[newvar.split('(')[0]+'_bar'] = alldic[vari+'_bar']
                    continue

    # Now change the system to include possible augmented variables
    endo_r = filter(lambda x: x[1] != [-1,0], spvdic2['endo']) 
    if endo_r:
        endo_r = [[x[0],[abs(x[1][0]+1),x[1][1]]] for x in endo_r ]
        # Create lags and forwards equations
        for vari in endo_r:
            for lag in range(abs(vari[1][0])):
                tind = (5-len(str(lag+2)))*'0'+str(lag+2)
                if lag == 0:
                    if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'(t)' not in list_tmp1:
                        list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'(t-1)')
                else:
                    tind1 = (5-len(str(lag+1)))*'0'+str(lag+1)
                    if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)' not in list_tmp1:
                        list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)')
            for lead in range(vari[1][1]):
                tind = (5-len(str(lead)))*'0'+str(lead)
                if lead == 0:
                    if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+vari[0][0] not in list_tmp1:
                        list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+vari[0][0])
                else:
                    tind1 = (5-len(str(lead-1)))*'0'+str(lead-1)
                    if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)' not in list_tmp1:
                        list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)')
    exo_r = filter(lambda x: x[1] != [0,1], spvdic2['exo'])
    if exo_r:
        exo_r = [[x[0],[abs(x[1][0]),x[1][1]-1]] for x in exo_r ]
        # Create lags and forwards equations
        for vari in exo_r:
            for lag in range(abs(vari[1][0])):
                tind = (5-len(str(lag+1)))*'0'+str(lag+1)
                if lag == 0:
                    if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'(t)' not in list_tmp1:
                        list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'(t)')
                else:
                    tind1 = (5-len(str(lag)))*'0'+str(lag)
                    if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)' not in list_tmp1:
                        list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)')
            if vari[1][1] > 0:
                for lead in range(vari[1][1]+1):
                    tind = (5-len(str(lead+1)))*'0'+str(lead+1)
                    if lead == 0:
                        if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'(t+1)' not in list_tmp1:
                            list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'(t+1)')
                    else:
                        tind1 = (5-len(str(lead)))*'0'+str(lead)
                        if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)' not in list_tmp1:
                            list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)')
    con_r = filter(lambda x: x[1] != [0,1], spvdic2['con']) 
    if con_r:
        con_r = [[x[0],[abs(x[1][0]),x[1][1]-1]] for x in con_r ]
        # Create lags and forwards equations
        for vari in con_r:
            for lag in range(abs(vari[1][0])):
                tind = (5-len(str(lag+1)))*'0'+str(lag+1)
                if lag == 0:
                    if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+'-'+vari[0][0].split('(')[0]+'(t)' not in list_tmp1:
                        list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+'-'+vari[0][0].split('(')[0]+'(t)')
                else:
                    tind1 = (5-len(str(lag)))*'0'+str(lag)
                    if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)' not in list_tmp1: 
                        list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)')
            if vari[1][1] > 0:
                for lead in range(int(vari[1][1])+1):
                    tind = (5-len(str(lead+1)))*'0'+str(lead+1)
                    if lead == 0:
                        if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'(t+1)' not in list_tmp1:
                            list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'(t+1)')
                    else:
                        tind1 = (5-len(str(lead)))*'0'+str(lead)
                        if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)' not in list_tmp1:
                            list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)')



    return self, (list_tmp1,list_tmp2)

def mkaug2(self, insys):
    # Determine the lengths of augmented vars

    list_tmp1 = deepcopy(insys)

    endosp = []
    for x in self.vardic['endo']['var']:
        endosp = endosp + [[x,[-1,0]]]
    exosp = []
    for x in self.vardic['exo']['var']:
        exosp = exosp + [[x,[0,1]]]
    consp = []
    for x in self.vardic['con']['var']:
        consp = consp + [[x,[0,1]]]

    alldic = {}
    alldic.update(self.paramdic)
    alldic.update(self.sstate)

    spvdic = {}
    spvdic['endo'] = endosp
    spvdic['exo'] = exosp
    spvdic['con'] = consp

    spvdic2 = deepcopy(spvdic)

    endoli = [x[0].split('(')[0].strip() for x in self.vardic['endo']['var']]
    exoli = [x[0].split('(')[0].strip() for x in self.vardic['exo']['var']]
    conli = [x[0].split('(')[0].strip() for x in self.vardic['con']['var']]


    patup = ('{-100,100}|None','all','{-100,100}')
    count = 0
    for i1,line in enumerate(list_tmp1):
        iterob = self.vreg(patup,line,True,'max')
        if iterob:
            iterob = list(iterob)
        else:
            continue
        iterob.reverse()
        for vma in iterob:
            vtype = vma[1][0]
            vartime = vma[2][2]
            vari = vma[2][1]
            pos = vma[3][0]
            poe = vma[3][1]
            if vtype == 'endo':
                indx = endoli.index(vari)
            elif vtype == 'con':
                indx = conli.index(vari)
            elif vtype == 'exo':
                indx = exoli.index(vari)
            else:
                continue

            # Check for endo
            if vtype == 'endo' and int(vartime) > spvdic['endo'][indx][1][1]:
                spvdic2['endo'][indx][1][1] = int(vartime)
                tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
                newvar = vari+'_F'+tind+'(t)'
                newname =  vari+'_F'+str(abs(int(vartime)))
                list_tmp1[i1] = list_tmp1[i1][:pos]+newvar.split('(')[0]+'(t)'+list_tmp1[i1][poe:]
                for i2 in range(int(vartime)):
                    tind = (5-len(str(abs(i2+1)-1)))*'0'+str(abs(i2+1)-1)
                    newvar = vari+'_F'+tind+'(t)'
                    newname =  vari+'_F'+str(abs(i2))
                    if [newvar,newname] not in self.vardic['con']['var']:
                        self.vardic['con']['var'].append([newvar,newname])
                        self.vardic['con']['mod'].append(self.vardic['endo']['mod'][indx])
                        self.audic['con']['var'].append([newvar,newname])
                        self.audic['con']['mod'].append(self.vardic['endo']['mod'][indx])
                        if self.sstate.has_key(vari+'_bar'):
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
                        else:
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
                    continue

            elif vtype == 'endo' and int(vartime) < spvdic['endo'][indx][1][0]:
                spvdic2['endo'][indx][1][0] = int(vartime)
                tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
                newvar = vari+'_B'+tind+'(t)'
                newname =  vari+'_B'+str(abs(int(vartime)))
                list_tmp1[i1] = list_tmp1[i1][:pos]+newvar.split('(')[0]+'(t-1)'+list_tmp1[i1][poe:]
                for i2 in range(1,abs(int(vartime))):
                    tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
                    newvar = vari+'_B'+tind+'(t)'
                    newname =  vari+'_B'+str(abs(i2+1))
                    if [newvar,newname] not in self.vardic['endo']['var']:
                        self.vardic['endo']['var'].append([newvar,newname])
                        self.vardic['endo']['mod'].append(self.vardic['endo']['mod'][indx])
                        self.audic['endo']['var'].append([newvar,newname])
                        self.audic['endo']['mod'].append(self.vardic['endo']['mod'][indx])
                        if self.sstate.has_key(vari+'_bar'):
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
                        else:
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
                    continue
            # Check for exo
            if vtype == 'exo' and int(vartime) > spvdic['exo'][indx][1][1]:
                spvdic2['exo'][indx][1][1] = int(vartime)
                tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
                newvar = vari+'_F'+tind+'(t)'
                newname =  vari+'_F'+str(abs(int(vartime)-1))
                list_tmp1[i1] = list_tmp1[i1][:pos]+newvar.split('(')[0]+'(t)'+list_tmp1[i1][poe:]
                for i2 in range(int(vartime)):
                    tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
                    newvar = vari+'_F'+tind+'(t)'
                    newname =  vari+'_F'+str(abs(i2+1))
                    if [newvar,newname] not in self.vardic['con']['var']:
                        self.vardic['con']['var'].append([newvar,newname])
                        self.vardic['con']['mod'].append(self.vardic['exo']['mod'][indx])
                        self.audic['con']['var'].append([newvar,newname])
                        self.audic['con']['mod'].append(self.vardic['exo']['mod'][indx])
                        if self.sstate.has_key(vari+'_bar'):
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
                        else:
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
                    continue
            elif vtype == 'exo' and int(vartime) < spvdic['exo'][indx][1][0]:
                spvdic2['exo'][indx][1][0] = int(vartime)
                tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
                newvar = vari+'_B'+tind+'(t)'
                newname =  vari+'_B'+str(abs(int(vartime)))
                list_tmp1[i1] = list_tmp1[i1][:pos]+newvar.split('(')[0]+'(t-1)'+list_tmp1[i1][poe:]
                for i2 in range(abs(int(vartime))):
                    tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
                    newvar = vari+'_B'+tind+'(t)'
                    newname =  vari+'_B'+str(abs(i2+1))
                    if [newvar,newname] not in self.vardic['endo']['var']:
                        self.vardic['endo']['var'].append([newvar,newname])
                        self.vardic['endo']['mod'].append(self.vardic['exo']['mod'][indx])
                        self.audic['endo']['var'].append([newvar,newname])
                        self.audic['endo']['mod'].append(self.vardic['exo']['mod'][indx])
                        if self.sstate.has_key(vari+'_bar'):
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
                        else:
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
                    continue
            # Check for con
            if vtype == 'con' and int(vartime) > spvdic['con'][indx][1][1]:
                spvdic2['con'][indx][1][1] = int(vartime)
                tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
                newvar = vari+'_F'+tind+'(t)'
                newname =  vari+'_F'+str(abs(int(vartime)-1))
                list_tmp1[i1] = list_tmp1[i1][:pos]+newvar.split('(')[0]+'(t)'+list_tmp1[i1][poe:]
                for i2 in range(int(vartime)):
                    tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
                    newvar = vari+'_F'+tind+'(t)'
                    newname =  vari+'_F'+str(abs(i2+1))
                    if [newvar,newname] not in self.vardic['con']['var']:
                        self.vardic['con']['var'].append([newvar,newname])
                        self.vardic['con']['mod'].append(self.vardic['con']['mod'][indx])
                        self.audic['con']['var'].append([newvar,newname])
                        self.audic['con']['mod'].append(self.vardic['con']['mod'][indx])
                        if self.sstate.has_key(vari+'_bar'):
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
                        else:
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
                    continue
            elif vtype == 'con' and int(vartime) < spvdic['con'][indx][1][0]:
                spvdic2['con'][indx][1][0] = int(vartime)
                tind = (5-len(str(abs(int(vartime)))))*'0'+str(abs(int(vartime)))
                newvar = vari+'_B'+tind+'(t)'
                newname =  vari+'_B'+str(abs(int(vartime)))
                list_tmp1[i1] = list_tmp1[i1][:pos]+newvar.split('(')[0]+'(t-1)'+list_tmp1[i1][poe:]
                for i2 in range(abs(int(vartime))):
                    tind = (5-len(str(abs(i2+2)-1)))*'0'+str(abs(i2+2)-1)
                    newvar = vari+'_B'+tind+'(t)'
                    newname =  vari+'_B'+str(abs(i2+1))
                    if [newvar,newname] not in self.vardic['endo']['var']:
                        self.vardic['endo']['var'].append([newvar,newname])
                        self.vardic['endo']['mod'].append(self.vardic['con']['mod'][indx])
                        self.audic['endo']['var'].append([newvar,newname])
                        self.audic['endo']['mod'].append(self.vardic['con']['mod'][indx])
                        if self.sstate.has_key(vari+'_bar'):
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.sstate[vari+'_bar']
                        else:
                            self.sstate[newvar.split('(')[0]+'_bar'] = self.paramdic[vari+'_bar']
                    continue


    # Now change the system to include possible augmented variables
    endo_r = filter(lambda x: x[1] != [-1,0], spvdic2['endo']) 
    if endo_r:
        endo_r = [[x[0],[abs(x[1][0]+1),x[1][1]]] for x in endo_r ]
        # Create lags and forwards equations
        for vari in endo_r:
            for lag in range(abs(vari[1][0])):
                tind = (5-len(str(lag+2)))*'0'+str(lag+2)
                if lag == 0:
                    if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'(t)' not in list_tmp1:
                        list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'(t-1)')
                else:
                    tind1 = (5-len(str(lag+1)))*'0'+str(lag+1)
                    if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)' not in list_tmp1:
                        list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)')
            for lead in range(vari[1][1]):
                tind = (5-len(str(lead)))*'0'+str(lead)
                if lead == 0:
                    if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+vari[0][0] not in list_tmp1:
                        list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+vari[0][0])
                else:
                    tind1 = (5-len(str(lead-1)))*'0'+str(lead-1)
                    if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)' not in list_tmp1:
                        list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)')
    exo_r = filter(lambda x: x[1] != [0,1], spvdic2['exo'])
    if exo_r:
        exo_r = [[x[0],[abs(x[1][0]),x[1][1]-1]] for x in exo_r ]
        # Create lags and forwards equations
        for vari in exo_r:
            for lag in range(abs(vari[1][0])):
                tind = (5-len(str(lag+1)))*'0'+str(lag+1)
                if lag == 0:
                    if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'(t)' not in list_tmp1:
                        list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'(t)')
                else:
                    tind1 = (5-len(str(lag)))*'0'+str(lag)
                    if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)' not in list_tmp1:
                        list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)')
            if vari[1][1] > 0:
                for lead in range(vari[1][1]+1):
                    tind = (5-len(str(lead+1)))*'0'+str(lead+1)
                    if lead == 0:
                        if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'(t+1)' not in list_tmp1:
                            list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'(t+1)')
                    else:
                        tind1 = (5-len(str(lead)))*'0'+str(lead)
                        if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)' not in list_tmp1:
                            list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)')
    con_r = filter(lambda x: x[1] != [0,1], spvdic2['con']) 
    if con_r:
        con_r = [[x[0],[abs(x[1][0]),x[1][1]-1]] for x in con_r ]
        # Create lags and forwards equations
        for vari in con_r:
            for lag in range(abs(vari[1][0])):
                tind = (5-len(str(lag+1)))*'0'+str(lag+1)
                if lag == 0:
                    if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+'-'+vari[0][0].split('(')[0]+'(t)' not in list_tmp1:
                        list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+'-'+vari[0][0].split('(')[0]+'(t)')
                else:
                    tind1 = (5-len(str(lag)))*'0'+str(lag)
                    if vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)' not in list_tmp1: 
                        list_tmp1.append(vari[0][0].split('(')[0]+'_B'+tind+'(t)'+' - '+vari[0][0].split('(')[0]+'_B'+tind1+'(t-1)')
            if vari[1][1] > 0:
                for lead in range(int(vari[1][1])+1):
                    tind = (5-len(str(lead+1)))*'0'+str(lead+1)
                    if lead == 0:
                        if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'(t+1)' not in list_tmp1:
                            list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'(t+1)')
                    else:
                        tind1 = (5-len(str(lead)))*'0'+str(lead)
                        if vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)' not in list_tmp1:
                            list_tmp1.append(vari[0][0].split('(')[0]+'_F'+tind+'(t)'+' - '+'E(t)|'+vari[0][0].split('(')[0]+'_F'+tind1+'(t+1)')

    return self, list_tmp1


def mk_subs_dic(self, secs):
    # Make the substitution system by joining lines and stripping
    list_tmp2 = []
    # This list is used to catch the @ALL instructions and store them
    list_tmp3 = []
    list_tmp3_ind = []
    # A list which collects all elements, to be handed to the template_paramdic for template parsing
    list_tmp4 = []
    i1=0
    linecounter=0
    for x in secs['vsfocs'][0]:
        if x.endswith(("...","\\")):
            linecounter += 1
        else:
            if linecounter == 0:
                ################# Catch @ALL expressions and store################
                if "@ALL" in x:
                    intor = []
                    intor1 = x.split('@ALL')[1].replace(';','')
                    intor.append('@ALL')
                    intor.append(intor1)
                    list_tmp4.append(intor)
                    raw_str = x.split(';')[0].split('@ALL')[1].split('{')[1].split('}')[0]
                    list_tmp3.append(raw_str.split(','))
                    list_tmp3_ind.append(i1)
                    continue
                ##################################################################
                line = x.replace(';','').split(']')[1].strip()
            elif linecounter > 0: # have a multiline equation
                str_tmp = ''
                for y in secs['vsfocs'][0][i1-linecounter:i1+1]:
                    str_tmp += y.replace('...','').replace(';','').replace('//','').strip()
                ################# Catch @ALL expressions and store################
                if "@ALL" in str_tmp:
                    intor = []
                    intor1 = x.split('@ALL')[1].replace(';','')
                    intor.append('@ALL')
                    intor.append(intor1)
                    list_tmp4.append(intor)
                    raw_str = str_tmp.split('@ALL')[1].split('{')[1].split('}')[0]
                    list_tmp3.append(raw_str.split(','))
                    list_tmp3_ind.append(i1)
                    continue
                ##################################################################
                if ']' in str_tmp:
                    line = str_tmp.split(']')[1].strip()
                else:
                    line = str_tmp.strip()
                linecounter = 0 
        i1 += 1
        splitline = line.split('=')
        list_tmp4.append([splitline[0].strip(), splitline[1].strip()])
        list_tmp2.append([splitline[0].strip(), splitline[1].strip()])
    self.allsubs_raw1 = deepcopy(list_tmp3)
    self.allsubs_index = deepcopy(list_tmp3_ind)
    self.nlsubs_raw1 = deepcopy(list_tmp2)
    self.nlsubsdic = dict(list_tmp2)
    if list_tmp4 != []: self.template_paramdic['subs_list'] = deepcopy(list_tmp4)
    else: self.template_paramdic['subs_list'] = False
    return self

def subs_in_subs_all(self):
    list_tmp3 = deepcopy(self.allsubs_raw1)
    list_tmp2 = deepcopy(self.nlsubs_raw1)
    # Replace substitutions inside substitutions, for both time-subscripted and steady state vars!
    mreg = re.compile('@DISCOUNT|@(E\(t.*?\)\|){0,1}.*?\(t.*?\)|@(E\(t.*?\)\|){0,1}.*?_bar')
    variables = [x[0] for x in list_tmp2] # do this once
    for i,x in enumerate(list_tmp3):
        expr = list_tmp3[i][0]
        while mreg.search(expr):
            ma = mreg.search(expr)
            indx = variables.index(ma.group())
            expr = list_tmp2[indx][1]
            list_tmp3[i][0] = expr
    self.allsubs_raw2 = deepcopy(list_tmp3)
    return self

def mk_all(self):
    list_tmp1 = deepcopy(self.allsubs_raw2)
    list_tmp2 = deepcopy(self.allsubs_raw1)
    list_tmp3 = deepcopy(self.nlsubs_raw1)
    subsli = [x[0] for x in list_tmp3]
    repldic_li = []
    do_ss = False
    for i1,elem in enumerate(list_tmp1):
        repldic = {}
        if 'SS' in elem: do_ss = True
        stem = list_tmp2[i1][0].split('(')[0]
        stem_time = '('+list_tmp2[i1][0].split('(')[1]
        # Do the simplest SS conversion
        if stem+'_bar' not in subsli: repldic[stem+'_bar'] = 'SS{'+stem+stem_time+'}'
        if '-' in elem[1]:
            chron = range(int(elem[1].split('-')[0].split('[')[1]),int(elem[1].split('-')[1].split(']')[0])+1)
        else:
            chron = elem[2].split('[')[1].split(']')[0].split(',')
            chron = [int(x) for x in chron]
        # Do current-period differentiation and SS jobs first
        patup = ('{-10,10}|None','endo|con|exo|other','{-10,10}')
        varli = self.vreg(patup,list_tmp1[i1][0],True,'max')
        for chrono in chron:
            if chrono == 0:
                for varo in varli:
                    if varo[1][0] == 'con' and int(varo[2][2]) == 0:
                        if stem+varo[2][1]+stem_time not in subsli:
                            repldic[stem+varo[2][1]+stem_time] = 'DIFF{'+stem+stem_time+','+varo[0]+'}'
                        if stem+varo[2][1]+'_bar' not in subsli:
                            repldic[stem+varo[2][1]+'_bar'] = 'SS{'+stem+varo[2][1]+stem_time+'}'
                    elif varo[1][0] == 'endo' and int(varo[2][2]) == -1:
                        if stem+varo[2][1]+stem_time not in subsli:
                            repldic[stem+varo[2][1]+stem_time] = 'DIFF{'+stem+stem_time+','+varo[0]+'}'
                        if stem+varo[2][1]+'_bar' not in subsli:
                            repldic[stem+varo[2][1]+'_bar'] = 'SS{'+stem+varo[2][1]+stem_time+'}'
            else:
                if stem_time == '(t)':
                    stem_time_new = '(t+'+str(chrono)+')'
                elif '+' in stem_time:
                    old_t = int(stem_time.split('+')[1].split(')')[0])
                    stem_time_new = '(t+'+str(old_t+chrono)+')'
                elif '-' in stem_time:
                    old_t = int('-'+stem_time.split('-')[1].split(')')[0])
                    stem_time_new = '(t+'+str(old_t+chrono)+')'
                for varo in varli:
                    if chrono > 0:
                        if varo[1][0] == 'con' and int(varo[2][2]) == 0:
                            if stem+stem_time_new not in subsli:
                                repldic[stem+stem_time_new] = 'FF_'+str(chrono)+'{'+stem+stem_time+'}'
                            if stem+varo[2][1]+stem_time_new not in subsli:
                                repldic[stem+varo[2][1]+stem_time_new] = 'DIFF{'+stem+stem_time_new+','+'FF_'+str(chrono)+'{'+varo[0]+'}'+'}'
                        elif varo[1][0] == 'endo' and int(varo[2][2]) == -1:
                            if stem+stem_time_new not in subsli:
                                repldic[stem+stem_time_new] = 'FF_'+str(chrono)+'{'+stem+stem_time+'}'
                            if stem+varo[2][1]+stem_time_new not in subsli:
                                repldic[stem+varo[2][1]+stem_time_new] = 'DIFF{'+stem+stem_time_new+','+'FF_'+str(chrono)+'{'+varo[0]+'}'+'}'
                    elif chrono < 0:
                        if varo[1][0] == 'con' and int(varo[2][2]) == 0:
                            if stem+stem_time_new not in subsli:
                                repldic[stem+stem_time_new] = 'BB_'+str(abs(chrono))+'{'+stem+stem_time+'}'
                            if stem+varo[2][1]+stem_time_new not in subsli:
                                repldic[stem+varo[2][1]+stem_time_new] = 'DIFF{'+stem+stem_time_new+','+'BB_'+str((chrono))+'{'+varo[0]+'}'+'}'
                        elif varo[1][0] == 'endo' and int(varo[2][2]) == -1:
                            if stem+stem_time_new not in subsli:
                                repldic[stem+stem_time_new] = 'BB_'+str(abs(chrono))+'{'+stem+stem_time+'}'
                            if stem+varo[2][1]+stem_time_new not in subsli:
                                repldic[stem+varo[2][1]+stem_time_new] = 'DIFF{'+stem+stem_time_new+','+'BB_'+str(abs(chrono))+'{'+varo[0]+'}'+'}'
        repldic_li.append(repldic)
    repldic_li.reverse()
    indexli = deepcopy(self.allsubs_index)
    indexli.reverse()
    for i1,indexo in enumerate(indexli):
        self.nlsubs_raw1 = self.nlsubs_raw1[:indexo]+[list(x) for x in repldic_li[i1].items()]+self.nlsubs_raw1[indexo:]
    del self.allsubs_index
    return self
                
            
            

def subs_in_subs(self):
    list_tmp2 = deepcopy(self.nlsubs_raw1)
    list_tmp3 = deepcopy(self.allsubs_raw1)
    # Replace substitutions inside substitutions, for both time-subscripted and steady state vars!
    mreg = re.compile('@DISCOUNT|@(E\(t.*?\)\|){0,1}.*?\(t.*?\)|@(E\(t.*?\)\|){0,1}.*?_bar')
    variables = [x[0] for x in list_tmp2] # do this once
    self.subs_vars = deepcopy(variables)
    for i,x in enumerate(list_tmp2):
        rhs_eq = list_tmp2[i][1]
        while mreg.search(rhs_eq):
            ma = mreg.search(rhs_eq)
            pos, poe = ma.span()
            indx = variables.index(ma.group())
            # Important: When substituting in, put the term between brackets!
            rhs_eq = rhs_eq[:pos]+'('+list_tmp2[indx][1]+')'+rhs_eq[poe:]
        # Don't do below anymore as the term above has been inserted between brackets, so okay that way
        '''
        # Finally get rid of possible `+-` or `-+` occurences
        while '+-' in rhs_eq or '-+' in rhs_eq or '++' in rhs_eq:
            rhs_eq = rhs_eq.replace('+-','-')
            rhs_eq = rhs_eq.replace('-+','-')
            rhs_eq = rhs_eq.replace('++','+')
        '''
        list_tmp2[i][1] = rhs_eq
    self.nlsubs_raw2 = deepcopy(list_tmp2)
    return self

def ff_chron_str(self,str1='',ff_int=1):
    '''
    This forwards all variables in an algebraic string expression by ff_int,
    but it leaves the timing of the expectations operator untouched.
    '''
    _mregv1b = '(?<!E)\(t(?P<oper>[\+|-]{0,1})(?P<counto>\d{0,2})\)'
    _mregv1c = '(?<=E)\(t(?P<oper>[\+|-]{0,1})(?P<counto>\d{0,2})\)(?=\|)'
    mregv1b = re.compile(_mregv1b)
    mregv1c = re.compile(_mregv1c)
    patup = ('{-10,10}|None','endo|con|exo|other','{-10,10}')
    var_li = list(self.vreg(patup,str1,True,'max'))
    var_li.reverse()
    for varo in var_li:
        varn = varo[0]
        ma = mregv1b.search(varn)
        starts = ma.start()
        ends = ma.end()
        matot = ma.group()
        oper = ma.group('oper')
        counto = ma.group('counto')
        varpos = varo[3]
        if oper == '':
            oper = '+'
            counto = int(ff_int)
        elif oper == '+':
            counto = int(counto)+ff_int
            counto = str(counto)
        elif oper == '-':
            counto = -int(counto)+ff_int
            counto = str(counto)
            if counto == '0':
                oper = ''
                counto = ''
            elif counto[0] == '-':
                oper == '-'
                counto = counto[1]
            else:
                oper = '+'
                counto = str(counto)
        # Split off expectations operator if present
        if '|' in varn:
            expos = varn.split('|')[0]
            rems = varn.split('|')[1]
        else:
            expos = ''
            rems = varn
        if expos == '':
            varn_new = rems.split('(')[0].replace('-','').replace('+','')+'(t'+str(oper)+str(counto)+')'
        else:
            varn_new = expos+'|'+rems.split('(')[0].replace('-','').replace('+','')+'(t'+str(oper)+str(counto)+')'
        # Add an expectations term if needed, but only if it is not already inside the original variable
        if oper == '+' and int(counto) > 0:
            if not mregv1c.search(varn_new): varn_new = 'E(t)|'+varn_new
        str1 = str1[:varpos[0]]+varn_new+str1[varpos[1]:]
    return str1

def bb_chron_str(self,str1='',bb_int=1):
    '''
    This backwards (lags) all variables in an algebraic string expression by bb_int,
    but it leaves the timing of the expectations operator untouched.
    '''
    _mregv1b = '(?<!E)\(t(?P<oper>[\+|-]{0,1})(?P<counto>\d{0,2})\)'
    mregv1b = re.compile(_mregv1b)    
    patup = ('{-10,10}|None','endo|con|exo|other','{-10,10}')
    var_li = list(self.vreg(patup,str1,True,'max'))
    var_li.reverse()
    for varo in var_li:
        varn = varo[0]
        ma = mregv1b.search(varn)
        starts = ma.start()
        ends = ma.end()
        matot = ma.group()
        oper = ma.group('oper')
        counto = ma.group('counto')
        varpos = varo[3]
        if oper == '':
            oper = '-'
            counto = int(bb_int)
        elif oper == '+':
            counto = int(counto)-bb_int
            counto = str(counto)
            if counto == '0':
                oper = ''
                counto = ''
            elif counto[0] == '-':
                oper == '-'
                counto = counto[1]
            else:
                oper = '+'
                counto = str(counto)
        elif oper == '-':
            counto = -int(counto)-bb_int
            counto = str(counto)[1]
        varn_new = varn.split('(')[0].replace('-','').replace('+','')+'(t'+str(oper)+str(counto)+')'
        # Compare this line with ff_chron_str, here I do not add any expectations term
        str1 = str1[:varpos[0]]+varn_new+str1[varpos[1]:]
    return str1

def ss_chron_str(self,str1=''):
    '''
    This turns all variables in an algebraic string expression into steady state equivalents.
    '''
    patup = ('{-10,10}|None','endo|con|exo|other','{-10,10}')
    try:
        var_li = list(self.vreg(patup,str1,True,'max'))
    except:
        print "Model parsing error, problem with processing this text string:\n"+"'"+str1+"'."
        sys.exit()

    var_li.reverse()
    for varo in var_li:
        varn = varo[0]
        varpos = varo[3]
        if '|' in varn:
            varn_new = varn.split('|')[1].split('(')[0]+'_bar'
        else:
            varn_new = varn.split('(')[0]+'_bar'
        str1 = str1[:varpos[0]]+varn_new+str1[varpos[1]:]
    return str1

def mk_steady(self):
    list_tmp2 = deepcopy(self.nlsubs_raw2)
    _mreg = 'SS\{.*?\}'
    mreg = re.compile(_mreg)
    for i1,elem in enumerate(list_tmp2):
        # If there is an unreplaced DIFF inside the SS, then skip for now...
        if list_tmp2[i1][1][:2] == 'SS' and 'DIFF' in list_tmp2[i1][1][2:] : continue
        while mreg.search(list_tmp2[i1][1]):
            ma = mreg.search(list_tmp2[i1][1])
            matn = ma.group()
            starts = ma.start()
            ends = ma.end()
            str_tmp = matn
            str_tmp = str_tmp.split('{')[1].split('}')[0]
            str_tmp = ss_chron_str(self,str1=str_tmp)
            list_tmp2[i1][1] = list_tmp2[i1][1][:starts]+str_tmp+list_tmp2[i1][1][ends:]
    self.nlsubs_raw2 = deepcopy(list_tmp2)  
    return self

def mk_steady2(self):
    list_tmp2 = deepcopy(self.foceqs2)
    patup = ('{-10,10}|None','endo|con|exo|other','{-10,10}')
    vreg = self.vreg
    for i1,elem in enumerate(list_tmp2):
        varli = list(vreg(patup,list_tmp2[i1],True,'max'))
        varli.reverse()
        for varo in varli:
            vpos = varo[3]
            list_tmp2[i1] = list_tmp2[i1][:vpos[0]]+varo[2][1]+'_bar'+list_tmp2[i1][vpos[1]:]
    self.foceqss = deepcopy(list_tmp2)  
    return self

def mk_timeshift(self):
    list_tmp2 = deepcopy(self.nlsubs_raw2)
    _mreg = '(?<=\*)FF[_](?P<fint>\d{1,2})\{.*?\}'
    _mreg2 = '(?<=\*)BB[_](?P<bint>\d{1,2})\{.*?\}'
    _mreg_b = '(?<!\*)FF[_](?P<fint>\d{1,2})\{.*?\}'
    _mreg2_b = '(?<!\*)BB[_](?P<bint>\d{1,2})\{.*?\}'
    _mreg3 = '\s*DIFF{.*?,.*?}\s*'
    mreg = re.compile(_mreg)
    mreg2 = re.compile(_mreg2)
    mregb= re.compile(_mreg_b)
    mreg2b = re.compile(_mreg2_b)
    mreg3 = re.compile(_mreg3)
    for i1,elem in enumerate(list_tmp2):
        # First search for forward shifting terms
        while mreg.search(list_tmp2[i1][1]):
            ma = mreg.search(list_tmp2[i1][1])
            matn = ma.group()
            # Skip if the is a part which needs to be differentiated
            if mreg3.search(matn): break            
            starts = ma.start()
            ends = ma.end()
            fint = int(ma.group('fint'))
            str_tmp = matn
            str_tmp = str_tmp.split('{')[1].split('}')[0]
            str_tmp = ff_chron_str(self,str1=str_tmp,ff_int=fint)
            list_tmp2[i1][1] = list_tmp2[i1][1][:starts]+str_tmp+list_tmp2[i1][1][ends:]
        # First search for forward shifting terms
        while mregb.search(list_tmp2[i1][1]):
            ma = mregb.search(list_tmp2[i1][1])
            matn = ma.group()
            # Skip if the is a part which needs to be differentiated
            if mreg3.search(matn): break            
            starts = ma.start()
            ends = ma.end()
            fint = int(ma.group('fint'))
            str_tmp = matn
            str_tmp = str_tmp.split('{')[1].split('}')[0]
            str_tmp = ff_chron_str(self,str1=str_tmp,ff_int=fint)
            list_tmp2[i1][1] = list_tmp2[i1][1][:starts]+str_tmp+list_tmp2[i1][1][ends:]
        # Now search for backward shifting terms
        while mreg2.search(list_tmp2[i1][1]):
            ma = mreg2.search(list_tmp2[i1][1])
            matn = ma.group()
            # Skip if the is a part which needs to be differentiated
            if mreg3.search(matn): break  
            starts = ma.start()
            ends = ma.end()
            bint = int(ma.group('bint'))
            str_tmp = matn
            str_tmp = str_tmp.split('{')[1].split('}')[0]
            str_tmp = bb_chron_str(self,str1=str_tmp,bb_int=bint)
            list_tmp2[i1][1] = list_tmp2[i1][1][:starts]+str_tmp+list_tmp2[i1][1][ends:]
        # Now search for backward shifting terms
        while mreg2b.search(list_tmp2[i1][1]):
            ma = mreg2b.search(list_tmp2[i1][1])
            matn = ma.group()
            # Skip if the is a part which needs to be differentiated
            if mreg3.search(matn): break  
            starts = ma.start()
            ends = ma.end()
            bint = int(ma.group('bint'))
            str_tmp = matn
            str_tmp = str_tmp.split('{')[1].split('}')[0]
            str_tmp = bb_chron_str(self,str1=str_tmp,bb_int=bint)
            list_tmp2[i1][1] = list_tmp2[i1][1][:starts]+str_tmp+list_tmp2[i1][1][ends:]
    self.nlsubs_raw2 = deepcopy(list_tmp2)  
    return self

def ext_differ_out(self):
    '''
    This function will check differentiation instructions and augment the original one
    for timeshifts and the discount factor
    '''
    list_tmp2 = deepcopy(self.nlsubs_raw2)
    # Also allow for differentiation in substitution list
    _mreg = '\s*DIFF{.*?,.*?}\s*'
    mreg = re.compile(_mreg)
    _mregv1 = '\w+\d*_bar'
    patup = ('{-10,10}|None','endo|con|exo|other','{-10,10}')
    for kk1,elem in enumerate(list_tmp2):
        # Skip any line in which the steady state needs to be take before DIFF
        if list_tmp2[kk1][1][:4] == 'DIFF' and 'SS{' in list_tmp2[kk1][1]: continue
        # Also skip cases such as @r_bar = DIFF{z_bar*k_bar**(rho)*l_bar**(1-rho),k_bar}
        if list_tmp2[kk1][1][:4] == 'DIFF' and not self.vreg(patup,list_tmp2[kk1][1],True,'max') and\
           "_bar" in list_tmp2[kk1][1].split(',')[1]: continue        
        if mreg.search(list_tmp2[kk1][1]):
            maoutli = mreg.finditer(list_tmp2[kk1][1])
            maoutli = [elem for elem in maoutli]
            maoutli.reverse()
            for maout in maoutli:
                expout = maout.group()
                starts_out = maout.start()
                ends_out = maout.end()
                evalstr = expout.replace('DIFF','')
                evalstr = evalstr.replace('{','')
                evalstr = evalstr.replace('}','')
                differo = evalstr.split(',')[1]
                evalstr = evalstr.split(',')[0]
                # When the variable differo is steady state then we can skip this test
                if '_bar' in differo: continue
                # Should be only one variable, so take the first element
                diff_li = list(self.vreg(patup,differo,True,'max'))[0]
                var_li = list(self.vreg(patup,evalstr,True,'max'))
                # Expectations timing and differos stem name
                diff_conf = diff_li[2][:2]
                # The differos time subscript
                diff_conf2 = int(diff_li[2][-1])
                # Check for case that we need to augment the expression to be differentiated
                modcase = False
                chrondiff = []
                for elem in var_li:
                    if elem[2][:2] == diff_conf and int(elem[2][-1]) != diff_conf2:
                        chrondiff.append((-1)*(diff_conf2+int(elem[2][-1])))
                        modcase = True
                if modcase:
                    if '@DISCOUNT' not in dict(self.nlsubs_raw2).keys():
                        print "ERROR: For this model you need to define the discount rate explicitly!"
                        print "You can do this by inserting the special reserved @DISCOUNT variable into the substitution section"
                        sys.exit()
                    else:
                        disco = dict(self.nlsubs_raw2)['@DISCOUNT']
                    maxfor = max(chrondiff)
                    minfor = min(chrondiff)
                    finalexp = 'DIFF{'+evalstr
                    for kk2 in range(minfor,maxfor+1):
                        if kk2 == 0: continue
                        elif kk2 > 0: finalexp += '+'+disco+'**'+'('+str(kk2)+')'+'*FF_'+str(kk2)+'{'+evalstr+'}'
                        elif kk2 < 0: finalexp += '+'+disco+'**'+'('+str(-kk2)+')'+'*BB_'+str(kk2)+'{'+evalstr+'}'
                    finalexp += ','+differo+'}'
                    list_tmp2[kk1][1] = list_tmp2[kk1][1][:starts_out]+finalexp+list_tmp2[kk1][1][ends_out+1:]
    self.nlsubs_raw2 = deepcopy(list_tmp2)
    return self
            


def differ_out(self):
    # Copy in raw nlsubs which has substitutions inside
    # substitutions already replaced by previous function call
    list_tmp2 = deepcopy(self.nlsubs_raw2)
    # Also allow for differentiation in substitution list
    _mreg = '\s*DIFF{.*?,.*?}\s*'
    mreg = re.compile(_mreg)
    _mregv1 = '\w+\d*_bar'
    _mregv1b = '(?<!E)\(t[\+|-]{0,1}\d{0,2}\)'
    mregv1b = re.compile(_mregv1b)
    _mregv1b2 = '(?<!E)__ll__t[\+|-]{0,1}\d{0,2}__rr__'
    mregv1b2 = re.compile(_mregv1b)
    _mregv1bb = '\(t[\+|-]\d{1,2}\)'    
    mregv1bb = re.compile(_mregv1bb)
    _mregv1bb2 = '__l__t[\+|-]\d{1,2}__r__'    
    mregv1bb2 = re.compile(_mregv1bb2)    
    _mregv1c = '\([p+|m+]\)'
    mregv1c = re.compile(_mregv1c)
    _mregv1d = 'E\(t[\+|-]{0,1}\d{0,2}\)\|'
    mregv1d = re.compile(_mregv1d)
    _mregv1e = 'E\[t[\+|-]{0,1}\d{0,2}]\|'
    mregv1e = re.compile(_mregv1e)
    _mreglog = 'LOG\(.*?\)'
    _mregexp = 'EXP\(.*?\)'
    mregle_ng = re.compile(_mregexp+'|'+_mreglog)
    _mreglog2 = 'LOG\[.*?]'
    _mregexp2 = 'EXP\[.*?]'
    mregle2_ng = re.compile(_mregexp2+'|'+_mreglog2)
    _mreglog = 'LOG\(.*\)'
    _mregexp = 'EXP\(.*\)'
    mregle_gg = re.compile(_mregexp+'|'+_mreglog)
    _mreglog2 = 'LOG\[.*]'
    _mregexp2 = 'EXP\[.*]'
    mregle2_gg = re.compile(_mregexp2+'|'+_mreglog2)
    _mregsh = '(?<!DI)FF[_]\d{1,2}|BB[_]\d{1,2}'
    mregsh = re.compile(_mregsh)
    patup = ('{-10,10}|None','endo|con|exo|other','{-10,10}')
    # Collect some necessary info for steady state substitutions
    var_bar = []
    varnames = deepcopy(self.vardic)
    for elem in ['con','endo','exo','other']:
        tmp_li = [x[0].split('(')[0]+'_bar' for x in varnames[elem]['var']]
        for elem2 in tmp_li:
            var_bar.append(elem2)
    for kk1,elem in enumerate(list_tmp2):
        # Skip any line in which
        if list_tmp2[kk1][1][:4] == 'DIFF' and 'SS{' in list_tmp2[kk1][1]: continue
        while mreg.search(list_tmp2[kk1][1]):
            maout = mreg.search(list_tmp2[kk1][1])
            expout = maout.group()
            starts_out = maout.end()
            ends_out = maout.start()
            evalstr = expout.replace('DIFF','')
            evalstr = evalstr.replace('{','')
            evalstr = evalstr.replace('}','')
            differo = evalstr.split(',')[1]
            # Skip if differo still has uncomputed timeshifter
            if mregsh.search(differo): break
            # Also skip if evalstr has uncomputed timeshifter
            if mregsh.search(evalstr): break
            evalstr = evalstr.split(',')[0]
            time_vars = False
            var_li = self.vreg(patup,evalstr,True,'max')
            if var_li: time_vars = True
            #########################################################
            # Do for steady state expressions, as is very easy, CASE1
            #########################################################
            if '_bar' in evalstr and not time_vars and '_bar' in differo:
                # Now substitute out exp and log in terms of sympycore expressions
                elog = re.compile('LOG\(')
                while elog.search(evalstr):
                    ma = elog.search(evalstr)
                    pos = ma.span()[0]
                    poe = ma.span()[1]
                    evalstr = evalstr[:pos]+'SP.log('+evalstr[poe:]
                eexp = re.compile('EXP\(')
                while eexp.search(evalstr):
                    ma = eexp.search(evalstr)
                    pos = ma.span()[0]
                    poe = ma.span()[1]
                    evalstr = evalstr[:pos]+'SP.exp('+evalstr[poe:]
                # Now populate scope with sympy symbols
                tmp_dic={}
                for elem in var_bar:
                    tmp_dic[elem] = SP.Symbol(elem)
                locals().update(tmp_dic)
                tmp_dic = {}
                for elem in self.paramdic.keys():
                    tmp_dic[elem] = SP.Symbol(elem)
                locals().update(tmp_dic)
                # Also expose any variables from the ssidic, just in case
                tmp_dic = {}
                for elem in self.ssidic.keys():
                    tmp_dic[elem] = SP.Symbol(elem)
                locals().update(tmp_dic)
                # Also expose the differo variable
                locals()[differo] = SP.Symbol(differo)
                # Population of scope done, now do calculation and continue in loop
                expr_bar = eval(evalstr)
                resstr = expr_bar.diff(locals()[str(differo)])
                list_tmp2[kk1][1] = str(resstr)
                continue
            ##############################################################################
            # Do for steady state expressions, as is very easy, CASE2, needs working varli
            ##############################################################################
            
            ############## Generate varli only here, as it was not needed above (it was, but only for testing) ############### 
            try:
                var_li = list(self.vreg(patup,evalstr,True,'max'))
                if var_li: time_vars = True
            except:
                print "Model parsing error, problem with processing this text string:\n"+"'"+list_tmp2[kk1][0]+" = "+list_tmp2[kk1][1]+"'."
                sys.exit()
            ###################################################################################################################

            if '_bar' not in evalstr and not time_vars and '_bar' in differo:
                # First replace all chronological variables in evalstr with _bar equivalents
                varbar_li = deepcopy(var_li)
                varbar_li.reverse()
                for varo in varbar_li:
                    vpos = varo[3]
                    evalstr = evalstr[:vpos[0]]+varo[2][1]+'_bar'+evalstr[vpos[1]:]
                # Now substitute out exp and log in terms of sympycore expressions
                elog = re.compile('LOG\(')
                while elog.search(evalstr):
                    ma = elog.search(evalstr)
                    pos = ma.span()[0]
                    poe = ma.span()[1]
                    evalstr = evalstr[:pos]+'SP.log('+evalstr[poe:]
                eexp = re.compile('EXP\(')
                while eexp.search(evalstr):
                    ma = eexp.search(evalstr)
                    pos = ma.span()[0]
                    poe = ma.span()[1]
                    evalstr = evalstr[:pos]+'SP.exp('+evalstr[poe:]
                # Now populate scope with sympy symbols
                tmp_dic={}
                for elem in var_bar:
                    tmp_dic[elem] = SP.Symbol(elem)
                locals().update(tmp_dic)
                tmp_dic = {}
                for elem in self.paramdic.keys():
                    tmp_dic[elem] = SP.Symbol(elem)
                locals().update(tmp_dic)
                # Also expose any variables from the ssidic, just in case
                tmp_dic = {}
                for elem in self.ssidic.keys():
                    tmp_dic[elem] = SP.Symbol(elem)
                locals().update(tmp_dic)
                # Also expose the differo variable
                locals()[differo] = SP.Symbol(differo)
                # Population of scope done, now do calculation and continue in loop
                expr_bar = eval(evalstr)
                resstr = expr_bar.diff(locals()[str(differo)])
                list_tmp2[kk1][1] = str(resstr)
                continue
            ###################################################
            # Steady State branch ends here
            ###################################################
            var_li2 = deepcopy(var_li)
            var_li2.reverse()
            ss_li = [x[2][1]+'_bar' for x in var_li2]
            # Replace t+1 or t-1 with something without the operators, in evalstr
            for jj1,elem in enumerate(var_li2):
                tmp_li = list(elem)
                # Replace expectations operators with something else inside evalstr
                if mregv1d.search(tmp_li[0]):
                    ma = mregv1d.search(tmp_li[0])
                    vname = ma.group()
                    ends = ma.end()
                    starts = ma.start()
                    vname = vname.replace('(','__ll__')
                    vname = vname.replace(')','__rr__')
                    vname = vname.replace('|','DELIM__')
                    tmp_li[0] = tmp_li[0][:starts]+vname+tmp_li[0][ends:]                      
                evalstr = evalstr[:elem[3][0]]+tmp_li[0].\
                replace('(','__l__').\
                replace(')','__r__').\
                replace('|','DELIM__')+\
                evalstr[elem[3][1]:]

            # Replace expectations operators with something else inside varli
            for i1,varo in enumerate(var_li):
                if mregv1d.search(varo[0]):
                    str1 = varo[0].split('|')[1]
                    str2 = varo[0].split('|')[0]
                    str2 = str2.replace('(','__ll__')
                    str2 = str2.replace(')','__rr__')
                    var_li[i1] = list(varo)
                    var_li[i1][0] = str2+'DELIM__'+str1
            # Replace the LOGs and EXPs with square brackets
            while mregle_ng.search(evalstr):
                ma = mregle_ng.search(evalstr)
                starts = ma.start()
                ends = ma.end()
                evalstr = evalstr[:ends-1] + ']' + evalstr[ends:]
                evalstr = evalstr[:starts+3] + '[' + evalstr[starts+4:]
            # Replace t+1 or t-1 with something without the operators, inside var_li
            for i1,varo in enumerate(var_li):
                var_li[i1] = list(varo)
                ma = mregv1bb.search(varo[0])
                if ma:
                    ends = ma.end()
                    starts = ma.start()
                    if '+' in varo[0]:
                        str1 = varo[0].split('(')[0]
                        counto = varo[0].split('+')[1][0]
                        var_li[i1][0] = str1+'__l__'+int(counto)*'p' + '__r__'
                    if '-' in varo[0]:
                        str1 = varo[0].split('(')[0]
                        counto = varo[0].split('-')[1][0]
                        var_li[i1][0] = str1+'__l__'+int(counto)*'m' + '__r__'
            while mregv1bb2.search(evalstr):
                ma = mregv1bb2.search(evalstr)
                varn = ma.group()
                ends = ma.end()
                starts = ma.start()
                if '+' in varn:
                    str1 = varn.split('__l__')[0]
                    counto = varn.split('+')[1][0]
                    varn_new = str1+'__l__'+int(counto)*'p' + '__r__'
                    evalstr = evalstr[:starts]+varn_new+evalstr[ends:]
                if '-' in varn:
                    str1 = varn.split('__l__')[0]
                    counto = varn.split('-')[1][0]
                    varn_new = str1+'__l__'+int(counto)*'m' + '__r__'
                    evalstr = evalstr[:starts]+varn_new+evalstr[ends:]
            # Replace left and right round brackets for variables and expose to sympycore
            for i1,varo in enumerate(var_li):
                str_tmp2 = var_li[i1][0]
                str_tmp2 = str_tmp2.replace('(','__l__')
                str_tmp2 = str_tmp2.replace(')','__r__')
                locals()[str_tmp2] = SP.Symbol(str_tmp2)
            for varo in self.paramdic.keys(): locals()[varo] = SP.Symbol(varo)
            # Replace left and right round brackets for the differo variable and expose to sympycore
            str_tmp3 = deepcopy(differo)
            # First, also check for (t+x),(t-x) or (t) in differo and replace with something readable by sympycore
            while mregv1b.search(str_tmp3):
                ma = mregv1b.search(str_tmp3)
                ends = ma.end()
                starts = ma.start()
                if '+' in str_tmp3:
                    counto = str_tmp3.split('+')[-1][0]
                    post = str_tmp3[starts:]
                    pre = str_tmp3[:starts]
                    # Also replace brackets
                    if ')' in pre: pre = pre.replace(')','__rr__')
                    if '(' in pre: pre = pre.replace('(','__ll__')
                    if ')' in post: post = post.replace(')','__r__')
                    if '(' in post: post = post.replace('(','__l__')                    
                    str_tmp3 = pre+'__l__'+int(counto)*'p' + '__r__' 
                elif '-' in str_tmp3:
                    counto = str_tmp3.split('-')[-1][0]
                    post = str_tmp3[starts:]
                    pre = str_tmp3[:starts]
                    # Also replace brackets
                    if ')' in pre: pre = pre.replace(')','__rr__')
                    if '(' in pre: pre = pre.replace('(','__ll__')
                    if ')' in post: post = post.replace(')','__r__')
                    if '(' in post: post = post.replace('(','__l__')                    
                    str_tmp3 = pre+'__l__'+int(counto)*'m' + '__r__'
                else:
                    post = str_tmp3[starts:]
                    pre = str_tmp3[:starts]
                    # Also replace brackets
                    if ')' in pre: pre = pre.replace(')','__rr__')
                    if '(' in pre: pre = pre.replace('(','__ll__')
                    if ')' in post: post = post.replace(')','__r__')
                    if '(' in post: post = post.replace('(','__l__')                    
                    str_tmp3 = pre+'__l__'+'t' + '__r__'                    
            differo_new = deepcopy(str_tmp3)
            # Also replace the delimiter when expectations are present in differo
            while '|' in differo_new: differo_new = differo_new.replace('|','DELIM__')
            locals()[differo_new] = SP.Symbol(differo_new)
            # Replace left and right brackets for the entire expression before calling diff
            var_li.reverse()
            for varo in var_li:
                varn_orig = varo[0]
                varn = varo[0].replace(')','__r__')
                varn = varn.replace('(','__l__')
                evalstr = evalstr.replace(varn_orig, varn)
            # Replace the LOGs and EXPs with round brackets again
            while mregle2_ng.search(evalstr):
                ma = mregle2_ng.search(evalstr)
                starts = ma.start()
                ends = ma.end()
                evalstr = evalstr[:ends-1] + ')' + evalstr[ends:]
                evalstr = evalstr[:starts+3] + '(' + evalstr[starts+4:]
            # Now substitute out exp and log in terms of sympycore expressions
            elog = re.compile('LOG\(')
            while elog.search(evalstr):
                ma = elog.search(evalstr)
                pos = ma.span()[0]
                poe = ma.span()[1]
                evalstr = evalstr[:pos]+'SP.log('+evalstr[poe:]
            eexp = re.compile('EXP\(')
            while eexp.search(evalstr):
                ma = eexp.search(evalstr)
                pos = ma.span()[0]
                poe = ma.span()[1]
                evalstr = evalstr[:pos]+'SP.exp('+evalstr[poe:]
            # Now evaluate and differentiate
            # Before evaluating expose the self.paramdic and self.ssidic to locals
            tmp_dic = {}
            for elem in self.paramdic.keys():
                tmp_dic[elem] = SP.Symbol(elem)
            locals().update(tmp_dic)
            # Also expose any variables from the ssidic, just in case
            if 'ssidic' in dir(self):
                tmp_dic = {}
                for elem in self.ssidic.keys():
                    tmp_dic[elem] = SP.Symbol(elem)
                locals().update(tmp_dic)
            # Also expose any variables from the ss_li, just in case
            tmp_dic = {}
            for elem in ss_li:
                tmp_dic[elem] = SP.Symbol(elem)
            locals().update(tmp_dic)
            expr = eval(evalstr)
            try:
                expr = eval(evalstr)
            except:
                print "Parse ERROR: Could not evaluate/differentiate expression: "
                print list_tmp2[kk1][0]+" = "+list_tmp2[kk1][1]
                sys.exit()
            resstr = expr.diff(locals()[str(differo_new)])
            resstr = str(resstr)
            # Now switch back to normal notation with brackets
            while '__l__' in resstr:
                resstr = resstr.replace('__l__','(')
            while '__r__' in resstr:
                resstr = resstr.replace('__r__',')')
            while '__ll__' in resstr:
                resstr = resstr.replace('__ll__','(')
            while '__rr__' in resstr:
                resstr = resstr.replace('__rr__',')')
            # Also replace the DELIM with |
            while 'DELIM__' in resstr: resstr = resstr.replace('DELIM__','|')
            # Now also switch back to normal t+1/t-x notation
            while mregv1c.search(resstr):
                ma = mregv1c.search(resstr)
                ends = ma.end()
                starts = ma.start()
                str3 = ma.group()
                lengo = len(str3[1:-1])
                if 'p' in str3: resstr = resstr[:starts]+'(t+'+str(lengo)+')'+resstr[ends:]
                elif 'm' in str3: resstr = resstr[:starts]+'(t-'+str(lengo)+')'+resstr[ends:]
            # Replace with the differentiate term, but also take brackets around it just in case    
            list_tmp2[kk1][1] = list_tmp2[kk1][1].replace(expout,'('+resstr+')')
            # Dont' do the below anymore as above I have put the differentiated term into brackets, so okay now that way
            '''
            # Finally get rid of possible `+(-` or `-(+` occurences
            while '+-' in list_tmp2[kk1][1] or '-+' in list_tmp2[kk1][1] or '++' in list_tmp2[kk1][1]:
                list_tmp2[kk1][1] = list_tmp2[kk1][1].replace('+(-','-(')
                list_tmp2[kk1][1] = list_tmp2[kk1][1].replace('-(+','-(')
                list_tmp2[kk1][1] = list_tmp2[kk1][1].replace('+(+','+(')
            '''
    self.nlsubs_raw2 = deepcopy(list_tmp2)
    self.nlsubs = deepcopy(dict(list_tmp2))
    self.nlsubs_list = deepcopy(list_tmp2)
    return self

def premknonlinsys(self,secs):
    '''
    This takes the raw lines of the FOCs and does some joining, splitting and stripping
    '''
    # Make the non-linear system by joining lines and stripping
    list_tmp1 = []
    i1 = 0
    linecounter = 0
    for x in secs['focs'][0]:
        if x.endswith(('...','\\')):
            linecounter += 1
        else:
            if linecounter == 0:
                if ']' in x: line = x.replace(';','').split(']')[1].strip()
                else: line = x.replace(';','').strip()
                line = line.split('=')[0].strip()
                list_tmp1.append(line)
            elif linecounter > 0:
                str_tmp = ''
                for y in secs['focs'][0][i1-linecounter:i1+1]:
                    if ']' in y: str_tmp += y.split(']')[1].replace('...','').replace('\\','').replace(';','').strip()
                    else: str_tmp += y.replace('...','').replace('\\','').replace(';','').strip()
                if ']' in x: line = str_tmp.split(']')[1].strip()
                else: line = str_tmp.strip()
                linecounter = 0 
                line = line.split('=')[0].strip()
                list_tmp1.append(line)
        i1 += 1
    self.foceqs = deepcopy(list_tmp1)
    return self

def premknonlinsys2(self, secs):
    list_tmp1 = deepcopy(self.foceqs)
    variables = deepcopy(self.subs_vars)
    i1 = 0    
    for x in list_tmp1:
        # substitute out in main nonlinear equation system
        list_tmp2 = deepcopy(self.nlsubs_list)
        mreg = re.compile('@(E\(t.*?\)\|){0,1}.*?\(t.*?\)')
        for i,x in enumerate(list_tmp1):
            rhs_eq = list_tmp1[i]
            while mreg.search(rhs_eq):
                ma = mreg.search(rhs_eq)
                subv = ma.group()
                pos, poe = ma.span()
                indx = variables.index(subv)
                rhs_eq = rhs_eq[:pos]+'('+list_tmp2[indx][1]+')'+\
                    rhs_eq[poe:]
            list_tmp1[i] = rhs_eq
        i1+=1
    self.foceqs2 = deepcopy(list_tmp1)
    return self

def mknonlinsys(self, secs):
    """
    Create Non-Linear FOC System
    """
    list_tmp1 = deepcopy(self.foceqs2)
    self, list_tmp1 = mkaug2(self, list_tmp1)

    if any([False if 'None' in x else True for x in secs['vsfocs'][0]]):
        variables = deepcopy(self.subs_vars)
        list_tmp2 = deepcopy(self.nlsubs_list)
        
        list_tmp3 = [x[1] for x in list_tmp2]
        self, outtup = mkaug1(self, list_tmp3,list_tmp1)

        list_tmp3 = outtup[1]
        list_tmp1 = outtup[0]
        for i1,x in enumerate(list_tmp3):
            list_tmp2[i1][1] = list_tmp3[i1]

        self.nlsubs = dict(list_tmp2)
        nlsubs2 = {}
        for x in [x[0] for x in self.vardic['other']['var']]:
            nlsubs2[x] = self.nlsubs['@'+x]
        self.nlsubs2 = nlsubs2

        # Create ordered nlsubsys
        if self.vardic['other']['var']:
            nlsubsys = []
            varother = self.vardic['other']['var']
            for vari in [x[0] for x in varother]:
                nlsubsys.append(nlsubs2[vari])
            self.nlsubsys = nlsubsys

    # Count the number of distinct forward expectations
    # AND the no of equations that containt them
    ffli = []
    count = 0
    patup = ('0','endo|con|exo','{1,10}')
    for line in list_tmp1[:-self.nexo]:
        outli = self.vreg(patup,line,True,'min')
        if outli:
            count = count + 1
            for elem in [x[0] for x in outli]:
                if elem not in ffli:
                    ffli.append(elem)
    self.ffen = count
    self.ffli = ffli
    self.fflin = len(ffli)

    # Now create simplified vdic with possibly new augmented vars
    vdic = dict([[x,self.vardic[x]['var']] for x in self.vardic.keys()])
    self.vdic = vdic

    # Finally, recomputed nendo, ncon, etc.
    self.nendo = len(self.vardic['endo']['var'])
    self.nexo = len(self.vardic['exo']['var'])
    self.ncon = len(self.vardic['con']['var'])
    self.nother = len(self.vardic['other']['var'])
    self.nstat = self.nendo+self.nexo
    self.nall = self.nstat+self.ncon

    self.nlsys_list = list_tmp1
    return self

def mkloglinsys2(inlist):
    """
    This function simply takes the left-hand side
    of each equation and takes it over to the
    right-hand side with a minus sign in front.
    Now we have equations of the form g(x)=0.
    Also takes care of fact that the term may
    already be negative. Also takes into account
    that some equations have already been written
    as g(x)=0!
    """
    _revar='(E\(t-{0,1}\+{0,1}\d*\)\|){0,1}[a-z]*(\(t-{0,1}\+{0,1}\d*\)){1,1}'
    re_var = re.compile(_revar)
    list_tmp1 = deepcopy(inlist)
    str_tmp2 = ''
    i1 = 0
    while i1 < len(list_tmp1):

        if list_tmp1[i1].split('=')[0].strip().strip(';') == '0':
            list_tmp1[i1] = list_tmp1[i1].split('=')[1].strip().strip(';')
            i1 = i1 + 1
            continue
        elif list_tmp1[i1].split('=')[1].strip().strip(';') == '0':
            list_tmp1[i1] = list_tmp1[i1].split('=')[0].strip().strip(';')
            i1 = i1 + 1
            continue

        list_tmp1[i1] = list_tmp1[i1].strip(';')[:]
        str_tmp2 = list_tmp1[i1].split('=')[0].strip()
        list_tmp1[i1] = list_tmp1[i1].split('=')[1].strip()

        while re_var.search(str_tmp2):
            ma1 = re_var.search(str_tmp2)
            pos1 = ma1.span()[0]
            poe1 = ma1.span()[1]
            # for coeff*var and coeff not begins with '-', but '+'
            if pos1 > 1 and str_tmp2[0] != '-' and str_tmp2[0] == '+':
                coeff = '-'+str_tmp2[1:pos1]
                vari = ma1.group(0)
                str_tmp2=str_tmp2[poe1:]
                if str_tmp2 == '':
                    str_tmp2 = '@'
            # for coeff*var and coeff not begins with '+', but '-'      
            elif pos1 > 1 and str_tmp2[0] != '+' and str_tmp2[0] == '-':
                coeff = '+'+str_tmp2[1:pos1]
                vari = ma1.group(0)
                str_tmp2=str_tmp2[poe1:]
                if str_tmp2 == '':
                    str_tmp2 = '@'
            # for coeff*var and coeff not begins with '+' and not begins with '-'
            elif pos1 > 1 and str_tmp2[0] != '+' and str_tmp2[0] != '-':
                coeff = '-'+str_tmp2[:pos1]
                vari = ma1.group(0)
                str_tmp2=str_tmp2[poe1:]
                if str_tmp2 == '':
                    str_tmp2 = '@'
            # for coeff == '-'
            elif pos1 == 1 and str_tmp2[0] == '-':
                coeff = '+'
                vari = ma1.group(0)
                str_tmp2=str_tmp2[poe1:]
            # for coeff == '+'
            elif pos1 == 1 and str_tmp2[0] == '+':
                coeff = '-'
                vari = ma1.group(0)
                str_tmp2=str_tmp2[poe1:]
            # for var*coeff
            elif pos1 == 0 and len(re_var.findall(str_tmp2)) > 1\
                 and str_tmp2[poe1] != '-' and str_tmp2[poe1] != '+':
                ma2 = re_var.search(str_tmp2[poe1:])
                pos2 = ma2.span()[0]+poe1
                coeff = '-'+str_tmp2[poe1:pos2]
                vari = ma1.group(0)
                str_tmp2=str_tmp2[pos2:]
            # for last bit
            elif pos1 == 0 and len(re_var.findall(str_tmp2)) == 1:
                coeff = '-'+str_tmp2[1:][poe1-1:]
                vari = ma1.group(0)
                str_tmp2=''
            # for coeff == ''
            elif pos1 == 0 and len(re_var.findall(str_tmp2)) > 1 and str_tmp2[poe1] == '-':
                coeff = '-'
                vari = ma1.group(0)
                str_tmp2=str_tmp2[poe1:]
            # for coeff == ''
            elif pos1 == 0 and len(re_var.findall(str_tmp2)) > 1 and str_tmp2[poe1] == '+':
                coeff = '-'
                vari = ma1.group(0)
                str_tmp2=str_tmp2[poe1:]
            # for last bit
            elif pos1 == 0 and len(re_var.findall(str_tmp2)) == 1:
                coeff = '-'
                vari = ma1.group(0)
                str_tmp2=''

            if coeff[-1] != '*' and coeff != '-':
                coeff = coeff+'*'
            list_tmp1[i1] = list_tmp1[i1]+coeff+vari
        str_tmp2 = ''   
        i1 = i1 + 1
    self.nlsys_list = deepcopy(list_tmp1)
    return list_tmp1


# Creates the log-linear system
def mkloglinsys1(secs):
    list_tmp = []
    i1=0
    counter=0
    for x in secs['modeq'][0]:
        if '...' in x:
            counter = counter + 1
        elif '...' not in x:
            if counter == 0:
                list_tmp.append(x.replace(';','').split(']')[1].strip())
            elif counter > 0:
                str_tmp = ''
                for y in secs['modeq'][0][i1-counter:i1+1]:
                    str_tmp = str_tmp + y.replace('...','').replace(';','').strip()
                list_tmp.append(str_tmp.split(']')[1].strip())
                counter = 0 
        i1=i1+1
    return list_tmp

# Create a Variance-Covariance matrix
def mksigmat(self, secs):
    str_tmp1 = []
    str_tmp1 = secs['vcvm'][0][0].split('[')[1].rstrip()
    if str_tmp1[-2:] == '];':
        str_tmp1=str_tmp1[:-2]
    if len(secs['vcvm']) > 1:
        for x in secs['vcvm'][0][1:]:
            if x[-2:] != '];':
                str_tmp1=str_tmp1+x.lstrip().rstrip()[:]
            elif x[-2:] == '];':
                str_tmp1=str_tmp1+x.lstrip().rstrip()[:-2]
    locals().update(self.paramdic)
    locals().update(self.sstate)
    list_tmp1 = str_tmp1.split(';')
    len1 = len(list_tmp1[0].split())
    mat1=mat.zeros((len1,len1))
    i1=0
    for x in list_tmp1:
        i2=0
        # Just another trap, only to be sure...
        if ']' in x: x = x.replace(']','')
        if '[' in x: x = x.replace('[','')
        for y in x.split():
            mat1[i1,i2] = float(eval(y))
            i2=i2+1
        i1=i1+1
    return mat1 

def mkeqtype(self):
    lsys = self.llsys_list
    tup1 = ('{-1,1}|None','iid|exo','{-1,1}')
    tup2 = ('{-1,1}|None','all','{-1,1}')
    tup3 = ('{-1,1}','all','{-1,1}')
    err_indx = []
    det_indx = []
    exp_indx = []
    for x,y in zip(lsys,range(len(lsys))):
        if self.vreg(('{-1,1}','all','{-1,1}'),x,False,'max'):
            exp_indx.append(y)
        elif self.vreg((None,'exo|iid','{-1,1}'),x,False,'max'):
            if len(self.vreg((None,'exo|iid','{-1,1}'),x,True,'max'))==\
               len(self.vreg(('{-1,1}|None','all','{-1,1}'),x,True,'max')):
                err_indx.append(y)
            else:
                det_indx.append(y)
        else:
            det_indx.append(y)

    self.eqindx = {}
    self.eqindx['err'] = err_indx
    self.eqindx['det'] = det_indx
    self.eqindx['exp'] = exp_indx
    return self

# Make symbolic system and numeric as well
def mksymsys(self):
    func = []
    subli = []
    patup = ('{-10,10}|None','all','{-10,10}')
    for x in self.sstate.keys()+self.paramdic.keys():
        locals()[x] = eval("SP.Symbol('"+x+"')")
    for y in self.llsys_list:
        str_tmp = y[:]
        vali = [x[0] for x in self.vreg(patup,y,True,'min')]
        vali2 = [[x,'sub'+str(u)] for x,u in zip(vali,range(0,len(vali),1))]
        vali3 = [(x[0],x[3]) for x in self.vreg(patup,str_tmp,True,'max')]
        vali3.reverse()
        valdic = dict(vali2)
        for x in vali3:
            str_tmp = str_tmp[:x[1][0]] + valdic[x[0]] + str_tmp[x[1][1]:]
        subli.append(valdic)
        for x in valdic.values():
            locals()[x] = eval("SP.Symbol('"+x+"')")
        func.append(eval(str_tmp))
    diffli = []
    i1 = 0
    for x in func:
        diffdic = {}
        for y in subli[i1].items():
            diffdic[y[0]] = eval('SP.diff(func[i1],'+y[1]+')')
        diffli.append(diffdic)
        i1=i1+1
    self.diffli1 = diffli
    diffli2 = []
    locals().update(self.sstate)
    locals().update(self.paramdic)
    for x in self.diffli1:
        tmpdic={}
        for y in x.keys():
#            tmpdic[y] = eval(x[y].tostr())
            tmpdic[y] = eval(str(x[y]))

        diffli2.append(tmpdic) 
    self.diffli2 = diffli2
    return self


#This function is needed in population stage 1, at the end
def mk_msstate_subs(self):
    """
    This function takes the manually defined numerical sstate conditions and then replaces
    any @terms with the corresponding substitutions
    """
    _mreg = '@[a-zA-Z]*_bar'
    mreg = re.compile(_mreg)
    tmp_list = deepcopy(self.ssys_list)
    sub_dic = deepcopy(self.nlsubs)
    for i1,x in enumerate(tmp_list):
        while mreg.search(tmp_list[i1]):
            ma = mreg.search(tmp_list[i1])
            str_tmp = ma.group()
            tmp_list[i1] = tmp_list[i1].replace(str_tmp,'('+sub_dic[str_tmp]+')')
    self.ssys_list = deepcopy(tmp_list)
    return self

#This function is needed in population stage 1, at the end
def mk_cfstate_subs(self):
    """
    This function takes the closed form defined numerical sstate conditions and then replaces
    any @terms with the corresponding substitutions
    """
    _mreg = '@[a-zA-Z]*_bar'
    mreg = re.compile(_mreg)
    tmp_list = deepcopy(self.manss_sys)
    sub_dic = deepcopy(self.nlsubs)
    for i1,x in enumerate(tmp_list):
        while mreg.search(tmp_list[i1]):
            ma = mreg.search(tmp_list[i1])
            str_tmp = ma.group()
            tmp_list[i1] = tmp_list[i1].replace(str_tmp,'('+sub_dic[str_tmp]+')')
    self.manss_sys = deepcopy(tmp_list)
    return self

# Extra population stage factored out, which is needed before steady state calculations
def populate_model_stage_one_a(self, secs):
    # Check and calculate the substitution list and dictionary
    if any([False if 'None' in x else True for x in secs['vsfocs'][0]]):
        # Create the raw nlsubs list 1
        self = mk_subs_dic(self, secs)       
        

# Extra population stage factored out, which is needed before steady state calculations
def populate_model_stage_one_b(self, secs):
    # Check and calculate the substitution list and dictionary
    if any([False if 'None' in x else True for x in secs['vsfocs'][0]]):
        # Substitute out in @ALL list for expressions, then do expressions
        self = subs_in_subs_all(self)
        self = mk_all(self)
        # Create the raw nlsubs list 2 by replacing substitutions inside substitutions
        self = subs_in_subs(self)
        # Apply steady state transformation where needed
        self = mk_steady(self)
        # Extend the differ out before doing the differ out
        self = ext_differ_out(self)
        # Apply timeshifts where needed
        self = mk_timeshift(self)
        # Make a first differentiation pass for DIFFs inside timeshifters
        self = differ_out(self)
        # Apply timeshifts where needed
        self = mk_timeshift(self)
        # Apply steady state transformation where needed
        self = mk_steady(self)
        # Extend the differ out before doing the differ out
        self = ext_differ_out(self)
        # Make second differentiation pass for remaining DIFFs
        self = differ_out(self)
        # Apply timeshifts where needed
        self = mk_timeshift(self)
    return self

def populate_model_stage_one_bb(self, secs):
    # Do substitutions inside the numerical steady state list
    # Check and do substitutions
    if all([False if 'None' in x else True for x in secs['vsfocs'][0]]) and\
       all([False if 'None' in x else True for x in secs['manualss'][0]]):
        # Do this only if USE_FOCS has not been used, otherwise ssys_list would be missing
        if '_internal_focs_used' not in dir(self): self = mk_msstate_subs(self)
    # Do substitutions inside the closed form steady state list
    # Check and do substitutions
    if all([False if 'None' in x else True for x in secs['vsfocs'][0]]) and\
       all([False if 'None' in x else True for x in secs['closedformss'][0]]):
        self = mk_cfstate_subs(self)
    # Prepare the nonlinear FOC system but only stripping the raw format
    if all([False if 'None' in x else True for x in secs['focs'][0]]):
        self = premknonlinsys(self,secs)
        # Save for template instantiation
        self.template_paramdic['focs_list'] = deepcopy(self.foceqs)
        if all([False if 'None' in x else True for x in secs['vsfocs'][0]]):
            # This takes the stripped system and does @ replacements
            self = premknonlinsys2(self,secs)
        else:
            self.foceqs2 = deepcopy(self.foceqs)
        # Also create a steady state version
        self = mk_steady2(self)
    else:
        # Save for template instantiation
        self.template_paramdic['focs_list'] = False
        
    return self


def populate_model_stage_two(self, secs):
    """
    2nd stage population of DSGE model.  Needs Steady State.
    """
    # Creates nonlinear foc system
    if any([False if 'None' in x else True for x in secs['focs'][0]]):
        self = mknonlinsys(self, secs)

    # Creates loglinear system
    if any([False if 'None' in x else True for x in secs['modeq'][0]]):
        llsys_list = mkloglinsys1(secs)
        self.llsys_list = mkloglinsys2(llsys_list)
        # Save for template instantiation
        self.template_paramdic['llsys_list'] = deepcopy(self.llsys_list)
        self = mksymsys(self) # creates symbolic and numerical system
        self = mkeqtype(self) # creates variance/covariance
    else:
        # Save for template instantiation
        self.template_paramdic['llsys_list'] = False

    # Creates variance/covariance
    if any([False if 'None' in x else True for x in secs['vcvm'][0]]) and\
       'sstate' in dir(self):
        self.sigma = mksigmat(self, secs)
        # Save for template instantiation
        self.template_paramdic['sigma'] = deepcopy(self.sigma)
    else:
        # Save for template instantiation
        self.template_paramdic['sigma'] = False

    return self
