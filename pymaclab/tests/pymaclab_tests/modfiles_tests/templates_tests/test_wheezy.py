import pymaclab as pm
import pymaclab.modfiles.models as models
import pymaclab.modfiles.templates.wheezy_template as template
import copy

# Instantiate model in the normal way
cee = pm.newMOD(models.testing.cee,mesg=True)
params = copy.deepcopy(cee.template_paramdic)

# Now use the template to do it again
modstr = template.render(params)
cee2 = pm.newMOD(modstr,mesg=True)

# Check if the two template_paramdics are identical
for keyo in cee.template_paramdic.keys():
    print "Now testing key: ",keyo
    if keyo == 'sigma':
        assert cee.template_paramdic[keyo].all() == cee2.template_paramdic[keyo].all()
    elif keyo == 'paramdic':
        for keyo2 in cee.template_paramdic[keyo].keys():
            assert round(cee.template_paramdic[keyo][keyo2],6) == round(cee2.template_paramdic[keyo][keyo2],6)
    else:
        assert cee.template_paramdic[keyo] == cee2.template_paramdic[keyo]


# Do for a couple of other models as well
modelli = []
for elem in dir(models.stable):
    if '__' not in elem: modelli.append(elem)
# Remove models which need external file or otherwise don't work
modelli.remove('rbc1_ext')
modelli.remove('jermann98')
modelli.remove('prog')

for modelo in modelli:
    print "Now testing for model: ",modelo
    exec("model1 = pm.newMOD(models.stable."+modelo+",mesg=True)")
    params = copy.deepcopy(model1.template_paramdic)
    modstr = template.render(params)
    model2 = pm.newMOD(modstr,mesg=True)
    # Check if the two template_paramdics are identical
    for keyo in model1.template_paramdic.keys():
        print "Now testing key: ",keyo
        if keyo == 'sigma':
            assert model1.template_paramdic[keyo].all() == model2.template_paramdic[keyo].all()
        elif keyo == 'paramdic':
            for keyo2 in model1.template_paramdic[keyo].keys():
                assert round(model1.template_paramdic[keyo][keyo2],6) == round(model2.template_paramdic[keyo][keyo2],6)
        else:
            assert model1.template_paramdic[keyo] == model2.template_paramdic[keyo]    


     
