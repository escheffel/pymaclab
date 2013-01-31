from copy import deepcopy
from pymaclab.dsge.translators import pml_to_dynarepp
from pymaclab.dsge.translators import dynarepp_to_pml
from pymaclab.dsge.translators import pml_to_pml


class Translators(object):
    def __init__(self,other=None):
        self.template_paramdic = deepcopy(other.template_paramdic)
    
    def pml_to_dynarepp(self,template_paramdic=None,fpath=None,focli=None):
        if fpath == None:
            if template_paramdic == None:
                return pml_to_dynarepp.translate(template_paramdic=self.template_paramdic,focli=focli)
            else:
                return pml_to_dynarepp.translate(template_paramdic=template_paramdic,focli=focli)
        else:
            if template_paramdic == None:
                pml_to_dynarepp.translate(template_paramdic=self.template_paramdic,fpath=fpath,focli=focli)
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