import pymaclab as pm
import pymaclab.modfiles.models as models

rbc = pm.newMOD(models.stable.rbc1_num,mesg=False,ncpus='auto')

def test_all():
    # Do for paramdic
    eta_key = 'eta'
    eta_old = 2.0
    eta_new = 5.0
    rho_key = 'rho'
    rho_old = 0.36
    rho_new = 0.35
    rbc.updaters_queued.paramdic[eta_key] = eta_new
    rbc.updaters_queued.paramdic[rho_key] = rho_new
    
    # Do for nlsubsdic
    U_key = '@U(t)'
    U_old = 'c(t)**(1-eta)/(1-eta)'
    U_new = 'c(t)**(1-eta*1.01)/(1-eta*1.01)'
    F_key = '@F(t)'
    F_old = 'z(t)*k(t-1)**rho'
    F_new = 'z(t)*k(t-1)**rho*1.01'
    rbc.updaters_queued.nlsubsdic[U_key] = U_new
    rbc.updaters_queued.nlsubsdic[F_key] = F_new
    
    # Do for vardic
    var_key = ['c(t)','consumption']
    indexo = rbc.vardic['con']['var'].index(var_key)
    var_old = 'bk'
    var_new = 'cf'
    rbc.updaters_queued.vardic['con']['mod'][indexo][1] = var_new
    
    
    # NOW process the queue
    rbc.updaters_queued.process_queue()
    
    
    # Did it work?
    # For paramdic
    assert rbc.updaters_queued.paramdic.wrapobj[eta_key] == eta_new
    assert rbc.updaters_queued.paramdic[eta_key] == eta_new
    assert rbc.paramdic[eta_key] == eta_new
    
    assert rbc.updaters_queued.paramdic.wrapobj[rho_key] == rho_new
    assert rbc.updaters_queued.paramdic[rho_key] == rho_new
    assert rbc.paramdic[rho_key] == rho_new
    
    # For nlsubsdic
    assert rbc.updaters_queued.nlsubsdic.wrapobj[U_key] == U_new
    assert rbc.updaters_queued.nlsubsdic[U_key] == U_new
    assert rbc.nlsubsdic[U_key] == U_new
    
    assert rbc.updaters_queued.nlsubsdic.wrapobj[F_key] == F_new
    assert rbc.updaters_queued.nlsubsdic[F_key] == F_new
    assert rbc.nlsubsdic[F_key] == F_new
    
    # For vardic
    assert rbc.updaters_queued.vardic.wrapobj['con']['mod'][indexo][1] == var_new
    assert rbc.updaters_queued.vardic['con']['mod'][indexo][1] == var_new
    assert rbc.vardic['con']['mod'][indexo][1] == var_new
    
    