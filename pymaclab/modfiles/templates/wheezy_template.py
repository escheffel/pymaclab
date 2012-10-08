from wheezy.template.engine import Engine as wheezyEngine
from wheezy.template.ext.core import CoreExtension as wheezyCoreExtension
from wheezy.template.loader import DictLoader as wheezyDictLoader


template = """\
@require(vardic,sigma,mod_desc,subs_list,focs_list,manss_sys,mod_name,mod_desc,llsys_list,paramdic,ssidic,ssys_list,use_focs)
%Model Description++++++++++++++++++++++++++++++++++
@if mod_desc:
@mod_desc;
@end
@if not mod_desc:
None
@end


%Model Information++++++++++++++++++++++++++++++++++
# Short model name
@if mod_name:
Name = @mod_name;
@end
# Short model description
@if mod_desc:
Description = @mod_desc;
@end



%Parameters+++++++++++++++++++++++++++++++++++++++++
@if paramdic:
@for item in [[x,y] for x,y in zip([[str(p[0]),str(p[1])] for p in paramdic.items()],[str(z+1) for z in range(len(paramdic.items()))])]:
[@item[1]]   @item[0][0] = @item[0][1];
@end
@end
@if not paramdic:
None
@end



%Variable Vectors++++++++++++++++++++++++++++++++++++
@if vardic:
@for x,y,z in zip(vardic['endo']['var'],vardic['endo']['mod'],[str(z+1) for z in range(len(vardic['endo']['var']))]):
@if len(y) == 2: 
[@z]   @x[0]:@x[1]{endo}[@y[0],@y[1]]
@else:
[@z]   @x[0]:@x[1]{endo}
@end
@end
@end
@if vardic:
@for x,y,z in zip(vardic['con']['var'],vardic['con']['mod'],[str(z+1) for z in range(len(vardic['con']['var']))]):
@if len(y) == 2: 
[@z]   @x[0]:@x[1]{con}[@y[0],@y[1]]
@else:
[@z]   @x[0]:@x[1]{con}
@end
@end
@end
@if vardic:
@for x,y,z in zip(vardic['exo']['var'],vardic['exo']['mod'],[str(z+1) for z in range(len(vardic['exo']['var']))]):
@if len(y) == 2: 
[@z]   @x[0]:@x[2]:@x[1]{exo}[@y[0],@y[1]]
@else:
[@z]   @x[0]:@x[2]:@x[1]{exo}
@end
@end
@end
@if vardic['other']['var'] != []:
@for x,y,z in zip([[v1,v2] for v1,v2 in zip(['@'+x[0] for x in vardic['other']['var']],[x[1] for x in vardic['other']['var']])],vardic['other']['mod'],[str(z+1) for z in range(len(vardic['other']['var']))]):
@if len(y) == 2: 
[@z]   @x[0]:@x[1] [@y[0],@y[1]]
@else:
[@z]   @x[0]:@x[1]
@end
@end
@end
@if not vardic:
None
@end


%Boundary Conditions+++++++++++++++++++++++++++++++++
None


%Variable Substitution Non-Linear System+++++++++++++
@if subs_list:
@for x,y in zip(subs_list,[str(z+1) for z in range(len(subs_list))]):
@if x[0] != '@ALL':
[@y]   @x[0] = @x[1];
@end
@if x[0] == '@ALL':
[@y]   @x[0]@x[1];
@end
@end
@end
@if not subs_list:
None
@end


%Non-Linear First-Order Conditions+++++++++++++++++++
@if focs_list:
@for x,y in zip(focs_list,[str(z+1) for z in range(len(focs_list))]):
[@y]   @x = 0;
@end
@end
@if not focs_list:
None
@end


%Steady States [Closed Form]++++++++++++++++++++++++++
@if manss_sys:
@for x,y in zip(manss_sys,[str(z+1) for z in range(len(manss_sys))]):
[@y]   @x
@end
@end
@if not manss_sys:
None
@end


%Steady State Non-Linear System [Manual]+++++++++++++
@if use_focs:
USE_FOCS=@str(use_focs);
@end
@if ssys_list:
@for x,y in zip(ssys_list,[str(z+1) for z in range(len(ssys_list))]):
[@y]   @x = 0;
@end
@end
@if not ssys_list and not use_focs:
None
@end

@if ssidic:
@for x,y in zip([[p1,str(p2)] for p1,p2 in ssidic.items()],[str(z+1) for z in range(len(ssidic.items()))]):
[@y]   @x[0] = @x[1];
@end
@end
@if not ssidic:
None
@end



%Log-Linearized Model Equations++++++++++++++++++++++
@if llsys_list:
@for x,y in zip(llsys_list,[str(z+1) for z in range(len(llsys_list))]):
[@y]   @x = 0;
@end
@end
@if not llsys_list:
None
@end


%Variance-Covariance Matrix++++++++++++++++++++++++++
@if type(sigma) != type(True):
Sigma = [@sigma.__str__().replace('[','').replace(']',';')[:-2] ];
@end
@else:
None
@end


%End Of Model File+++++++++++++++++++++++++++++++++++

"""

engine = wheezyEngine(
    loader=wheezyDictLoader({'wheezy': template}),
    extensions=[wheezyCoreExtension()]
)
wheezy_template = engine.get_template('wheezy')