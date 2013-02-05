from copy import deepcopy
from pymaclab.dsge.translators import pml_to_dynarepp
from pymaclab.dsge.translators import dynarepp_to_pml
from pymaclab.dsge.translators import pml_to_pml
from pymaclab.dsge.parsers._dsgeparser import ff_chron_str, bb_chron_str


class Translators(object):
    def __init__(self,other=None):
        self.template_paramdic = deepcopy(other.template_paramdic)
        self._other = other
        
    
    def pml_to_dynarepp(self,template_paramdic=None,fpath=None,focli=None):
        # Need to do some work to make focs_li dynare-conformable!
        focli = self.template_paramdic['focs_dynare']
        other = self._other
        vreg = other.vreg
        patup = ('{-10,10}|None','endo|con|exo|iid|other','{-10,10}')
        compset = set(['endo','con'])
        for i1,lino in enumerate(focli):
            varli = set([x[1][0] for x in vreg(patup,focli[i1],True,'max')])
            if varli.intersection(compset) != set([]) and 'exo' in varli:
                focli[i1] = ff_chron_str(other,str1=focli[i1],ff_int=1,vtype='exo')
            else:
                focli[i1] = bb_chron_str(other,str1=focli[i1],bb_int=1,vtype='iid')
        template_paramdic = deepcopy(self.template_paramdic)
        template_paramdic['focs_dynare'] = focli
        
        if fpath == None:
            return pml_to_dynarepp.translate(template_paramdic=template_paramdic,focli=focli)
        else:
            pml_to_dynarepp.translate(template_paramdic=template_paramdic,fpath=fpath,focli=focli)
    
    def dynarepp_to_pml(self,template_paramdic=None,fpath=None,focli=None):
        if fpath == None:
            if template_paramdic == None:
                return dynarepp_to_pml.translate(template_paramdic=self.template_paramdic,focli=focli)
            else:
                return dynarepp_to_pml.translate(template_paramdic=template_paramdic,focli=focli)
        else:
            if template_paramdic == None:
                dynarepp_to_pml.translate(template_paramdic=self.template_paramdic,fpath=fpath,focli=focli)
            else:
                dynarepp_to_pml.translate(template_paramdic=template_paramdic,fpath=fpath,focli=focli)

    def pml_to_pml(self,template_paramdic=None,fpath=None):
        if fpath == None:
            if template_paramdic == None:
                return pml_to_pml.translate(template_paramdic=self.template_paramdic)
            else:
                return pml_to_pml.translate(template_paramdic=template_paramdic)
        else:
            if template_paramdic == None:
                pml_to_pml.translate(template_paramdic=self.template_paramdic,fpath=fpath)
            else:
                pml_to_pml.translate(template_paramdic=template_paramdic,fpath=fpath)