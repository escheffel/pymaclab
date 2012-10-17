from mako.template import Template


template = """\

var ${[x[0].split('(')[0] for x in vardic['endo']['var']+vardic['con']['var']+vardic['exo']['var']].__str__().replace('[','').replace(']','').replace("'","")};
varexo ${[x[2].split('(')[0] for x in vardic['exo']['var']].__str__().replace('[','').replace(']','').replace("'","")};

parameters ${[x for x in paramdic.keys() if '_bar' not in x].__str__().replace('[','').replace(']','').replace("'","")};

% for keyo in [x for x in paramdic.keys() if '_bar' not in x]:
${keyo} = ${paramdic[keyo]};
% endfor


model;
% for item in focs_list2:
0=${item.replace('E(t)|','').replace('LOG','log').replace('EXP','exp').replace('(t)','').replace('(t+1)','(+1)').replace('(t-1)','(-1)').replace('**','^')};
% endfor
end;

<%
initvli = []
allvar = [x[0].split('(')[0] for x in vardic['endo']['var']+vardic['exo']['var']+vardic['con']['var']]
%>
initval;
% for item in ssili:
% if item[0] in [x[0].split('(')[0]+'_bar' for x in vardic['con']['var']] or item[0] in [x[0].split('(')[0]+'_bar' for x in vardic['endo']['var']] or item[0] in [x[0].split('(')[0]+'_bar' for x in vardic['exo']['var']]:
${item[0].replace('_bar','')} = ${item[1]}; <% initvli.append(item[0].replace('_bar','')) %>
% endif
% endfor
% for keyo in [x for x in paramdic.keys() if '_bar' in x]:
% if keyo in [x[0].split('(')[0]+'_bar' for x in vardic['con']['var']] or keyo in [x[0].split('(')[0]+'_bar' for x in vardic['endo']['var']] or keyo in [x[0].split('(')[0]+'_bar' for x in vardic['exo']['var']]:
${keyo.replace('_bar','')} = ${paramdic[keyo]}; <% initvli.append(item[0].replace('_bar','')) %>
% endif
% endfor
% for varo in [x for x in [y for y in allvar if y not in initvli]]:
${varo} = 1.0;
% endfor
end;


order = 1;

% for i1,row in enumerate(range(sigma.shape[0])):
% if i1 == 0 and sigma.shape[0] == 1:
vcov = [${sigma[i1].__str__().replace('[','').replace(']','').replace("'","")} ];
% elif i1 == 0 and sigma.shape[0] != 1:
vcov = [${sigma[i1].__str__().replace('[','').replace(']','').replace("'","")};
% elif i1 != 0 and i1 != sigma.shape[0]-1:
${sigma[i1].__str__().replace('[','').replace(']','').replace("'","")};
% elif i1 != 0 and i1 == sigma.shape[0]-1:
${sigma[i1].__str__().replace('[','').replace(']','').replace("'","")} ];
% endif
% endfor


"""

mako_dynare_template = Template(template)