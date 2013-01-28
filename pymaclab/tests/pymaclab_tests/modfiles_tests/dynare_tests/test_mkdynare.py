import pymaclab as pm
import pymaclab.modfiles.models as models
import pymaclab.modfiles.templates.wheezy_template as template
import copy


def test_mkdynare_basic():

    # Test mk_dynare method on model
    rbc = pm.newMOD(models.stable.rbc1_focs,mesg=True)
    rbc.modsolvers.dynarepp.solve()
    
    # Test mk_dynare method on model
    rbc = pm.newMOD(models.stable.rbc1_cf,mesg=True)
    rbc.modsolvers.dynarepp.solve()
    
    # Test mk_dynare method on model
    rbc = pm.newMOD(models.stable.rbc1_num,mesg=True)
    rbc.modsolvers.dynarepp.solve()
    
    # Test mk_dynare method on model
    rbc = pm.newMOD(models.stable.rbc1_res,mesg=True)
    rbc.modsolvers.dynarepp.solve()
    
    # Test mk_dynare method on model
    rbc = pm.newMOD(models.stable.rbc1_sug,mesg=True)
    rbc.modsolvers.dynarepp.solve()
    
def test_mkdynare_other():

    '''
    print "Now testing dynare translation of Monika Merz's model"
    # Test mk_dynare method on model
    rbc = pm.newMOD(models.stable.merz,mesg=True)
    rbc.mk_dynare()
    '''
    '''
    print "Now testing dynare translation of RBC2 model"
    # Test mk_dynare method on model
    rbc = pm.newMOD(models.stable.rbc2,mesg=True)
    rbc.mk_dynare()
    '''
    
    print "Now testing dynare translation of Grohe and Uribe's RBC model"
    # Test mk_dynare method on model
    rbc = pm.newMOD(models.stable.grohurib03,mesg=True)
    rbc.modsolvers.dynarepp.solve()
    
    print "Now testing dynare translation of Chris Sim's RBC model"
    # Test mk_dynare method on model
    rbc = pm.newMOD(models.stable.sims,mesg=True)
    rbc.modsolvers.dynarepp.solve()
    



     
