"""
Functions to parse the MODparser and attach results to the DSGE Model.
"""
import re
from copy import deepcopy
import numpy.matlib as mat
import sympycore as SP

#TODO: why have these as nested functions?
def populate_model_stage_one(self, secs):
    """
    1st stage population of DSGE model.  Does not need Steady State.
    """

#    _nreg = '^\s*None\s*$'
#    nreg = re.compile(_nreg)

    # Get all information from the varvec section
    if any([False if 'None' in x else True for x in secs['varvec'][0]]):
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

    # Extract the model name and description
    if any([False if 'None' in x else True for x in secs['info'][0]]):
       for x in secs['info'][0]:
            if 'Name' in x:
                self.mod_name = x.split('=')[1].replace(';','').strip()
            if 'Desc' in x:
                self.mod_desc = x.split('=')[1].replace(';','').strip()

    # Extract parameters into dictionary
    if any([False if 'None' in x else True for x in secs['params'][0]]):
        param = {}
        # need to do this so users don't have to worry about integer division
        # but preserves integer division for sympycore stuff
        safe_div = """
from __future__ import division
"""
        for x in secs['params'][0]:
            list_tmp = x.split(';')
            list_tmp = list_tmp[0].split('=')[:]
            str_tmp1 = list_tmp[0].strip()
            str_tmp2 = list_tmp[1].strip()
            #NOTE: this *should* be safe, but users should know what's
            # in the .mod file
            exec(safe_div + "param['"+str_tmp1+"']=" + str_tmp2, {}, locals())
            locals()[str_tmp1] = param[str_tmp1]
        self.paramdic = param

    # Collect information on manual (closed-form) steady state
    if any([False if 'None' in x else True for x in secs['closedformss'][0]]):
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
                    list_tmp1.append(x)
                elif counter > 0:
                    str_tmp = ''
                    for y in mansys[i1-counter:i1+1]:
                        str_tmp = str_tmp + y.replace('...','')
                    list_tmp1.append(str_tmp)
                    counter = 0 
            i1=i1+1
        self.manss_sys = list_tmp1

    # Collect info on numerical steady state
    if any([False if 'None' in x else True for x in secs['manualss'][0]]):
        _mreg = '[a-zA-Z]*_bar\s*=\s*[0-9]*\.[0-9]*'
        _mreg2 = '[a-zA-Z]*\s*=\s*[0-9]*\.[0-9]*'
        mreg = re.compile(_mreg+'|'+_mreg2)
        indx = []
        ssidic={}
        list_tmp = []
        i1=0
        counter=0
        for x in secs['manualss'][0]:
            if mreg.search(x):
                ma = mreg.search(x)
                str1 = ma.group().split('=')[0].strip()
                str2 = ma.group().split('=')[1].strip()
                ssidic[str1] = eval(str2)
                indx.append(i1)
            elif not mreg.search(x) and '...' in x:
                counter = counter + 1
            elif not mreg.search(x) and '...' not in x:
                if counter == 0:
                    list_tmp.append(x.replace(';','').split(']')[1].split('=')[0].strip())
                elif counter > 0:
                    str_tmp = ''
                    for y in secs['manualss'][0][i1-counter:i1+1]:
                        str_tmp = str_tmp + y.replace('...','').strip()
                    list_tmp.append(str_tmp.split(']')[1].split('=')[0].replace(';','').strip())
                    counter = 0 
            i1=i1+1
        self.ssidic = ssidic
        self.ssys_list = list_tmp

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


def mknonlinsys(self, secs):
    """
    Create Non-Linear FOC System
    """
    # Make the non-linear system by joining lines and stripping
    list_tmp1 = []
    i1 = 0
    linecounter = 0
    for x in secs['focs'][0]:
        if x.endswith(('...','\\')):
            linecounter += 1
        else:
            if linecounter == 0:
                line = x.replace(';','').split(']')[1].strip()
                line = line.split('=')[0].strip()
                list_tmp1.append(line)
            elif linecounter > 0:
                str_tmp = ''
                for y in secs['focs'][0][i1-linecounter:i1+1]:
                    str_tmp += y.replace('...','').replace('\\',
                                    '').replace(';','').strip()
                line = str_tmp.split(']')[1].strip()
                linecounter = 0 
                line = line.split('=')[0].strip()
                list_tmp1.append(line)
        i1 += 1

    # Check and do substitutions
    if any([False if 'None' in x else True for x in secs['vsfocs'][0]]):
        # Make the substitution system by joining lines and stripping
        list_tmp2 = []
        i1=0
        linecounter=0
        for x in secs['vsfocs'][0]:
            if x.endswith(("...","\\")):
                linecounter += 1
            else:
                if linecounter == 0:
                    line = x.replace(';','').split(']')[1].strip()
                elif linecounter > 0: # have a multiline equation
                    str_tmp = ''
                    for y in secs['vsfocs'][0][i1-linecounter:i1+1]:
                        str_tmp += y.replace('...','').replace(';',
                                        '').replace('//','').strip()
                    line = str_tmp.split(']')[1].strip()
                    linecounter = 0 
            i1 += 1
            splitline = line.split('=')
            list_tmp2.append([splitline[0].strip(), splitline[1].strip()])

        # Replace substitutions inside substitutions
        mreg = re.compile('@(E\(t.*?\)\|){0,1}.*?\(t.*?\)')
        variables = [x[0] for x in list_tmp2] # do this once
        for i,x in enumerate(list_tmp2):
            rhs_eq = list_tmp2[i][1]
            while mreg.search(rhs_eq):
                ma = mreg.search(rhs_eq)
                pos, poe = ma.span()
                indx = variables.index(ma.group())
                rhs_eq = rhs_eq[:pos]+'('+list_tmp2[indx][1]+')'+\
                                        rhs_eq[poe:]
            list_tmp2[i][1] = rhs_eq


        # substitute out in main nonlinear equation system
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

    self, list_tmp1 = mkaug2(self, list_tmp1)
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
        self = mksymsys(self) # creates symbolic and numerical system
        self = mkeqtype(self) # creates variance/covariance

    # Creates variance/covariance
    if any([False if 'None' in x else True for x in secs['vcvm'][0]]) and\
            'sstate' in dir(self):
        self.sigma = mksigmat(self, secs)

    return self




