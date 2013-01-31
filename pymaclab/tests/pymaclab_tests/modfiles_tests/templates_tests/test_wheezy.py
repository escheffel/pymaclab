import pymaclab as pm
import pymaclab.modfiles.models as models
import copy


def test_cee():

    # Instantiate model in the normal way
    cee = pm.newMOD(models.testing.cee,mesg=True)
    params = copy.deepcopy(cee.template_paramdic)

    # Now use the template to do it again
    modstr = cee.translators.pml_to_pml(template_paramdic=params)
    cee2 = pm.newMOD(modstr,mesg=True)

    # Check if the two template_paramdics are identical
    for keyo in cee.template_paramdic.keys():
        print "Now testing key: ",keyo
        print cee.template_paramdic[keyo]
        print '--------------------------'
        print cee2.template_paramdic[keyo]
        if keyo == 'sigma':
            assert cee.template_paramdic[keyo].all() == cee2.template_paramdic[keyo].all()
        elif keyo == 'paramdic':
            for keyo2 in cee.template_paramdic[keyo].keys():
                assert round(cee.template_paramdic[keyo][keyo2],6) == round(cee2.template_paramdic[keyo][keyo2],6)
        elif keyo == 'ssili':
            for i1,elem in enumerate([x[0] for x in cee.template_paramdic[keyo]]):
                assert round(float(cee.template_paramdic[keyo][i1][1]),6) == round(float(cee2.template_paramdic[keyo][i1][1]),6)
        elif keyo == 'ssidic':
            for keyo2 in cee.template_paramdic[keyo].keys():
                assert round(cee.template_paramdic[keyo][keyo2],6) == round(cee2.template_paramdic[keyo][keyo2],6)
        else:
            assert cee.template_paramdic[keyo] == cee2.template_paramdic[keyo]


def test_others():

    # Do for a couple of other models as well
    modelli = []
    for elem in dir(models.stable):
        if '__' not in elem: modelli.append(elem)
    # Remove models which need external file or otherwise don't work
    modelli.remove('rbc1_ext')
    modelli.remove('rbc1_extss')
    modelli.remove('jermann98_ext')
    modelli.remove('jermann98')
    modelli.remove('grohurib03')
    modelli.remove('merz')

    for modelo in modelli:
        print "Now testing for model: ",modelo
        exec(modelo+" = pm.newMOD(models.stable."+modelo+",mesg=True)")
        params = copy.deepcopy(eval(modelo+".template_paramdic"))
        exec('modstr = '+modelo+'.translators.pml_to_pml(template_paramdic=params)')
        exec(modelo+'_alt'+" = pm.newMOD(modstr,mesg=True)")
        # Check if the two template_paramdics are identical
        for keyo in eval(modelo+".template_paramdic.keys()"):
            print "Now testing key: ",keyo
            if keyo == 'sigma':
                assert eval(modelo+".template_paramdic[keyo].all()") == eval(modelo+"_alt"+".template_paramdic[keyo].all()")
            elif keyo == 'paramdic':
                for keyo2 in eval(modelo+".template_paramdic[keyo].keys()"):
                    assert round(eval(modelo+".template_paramdic[keyo][keyo2]"),6) == round(eval(modelo+"_alt"+".template_paramdic[keyo][keyo2]"),6)
            elif keyo == 'ssidic' and eval(modelo+".template_paramdic[keyo]") != False:
                #assert eval(modelo+".template_paramdic[keyo].keys()") == eval(modelo+"_alt"+".template_paramdic[keyo].keys()")
                for keyo2 in eval(modelo+".template_paramdic[keyo].keys()"):
                    try:
                        float(eval(modelo+".template_paramdic[keyo][keyo2]"))
                        assert round(eval(modelo+".template_paramdic[keyo][keyo2]"),6) == round(eval(modelo+"_alt"+".template_paramdic[keyo][keyo2]"),6)
                    except:
                        eval(modelo+".template_paramdic[keyo][keyo2]") == eval(modelo+"_alt"+".template_paramdic[keyo][keyo2]")
            elif keyo == 'paramdic':
                #assert eval(modelo+".template_paramdic[keyo].keys()") == eval(modelo+"_alt"+".template_paramdic[keyo].keys()")
                for keyo2 in eval(modelo+".template_paramdic[keyo].keys()"):
                    try:
                        float(eval(modelo+".template_paramdic[keyo][keyo2]"))
                        assert round(eval(modelo+".template_paramdic[keyo][keyo2]"),6) == round(eval(modelo+"_alt"+".template_paramdic[keyo][keyo2]"),6)
                    except:
                        eval(modelo+".template_paramdic[keyo][keyo2]") == eval(modelo+"_alt"+".template_paramdic[keyo][keyo2]")



