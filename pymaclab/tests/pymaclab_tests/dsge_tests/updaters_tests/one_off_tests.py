import pymaclab as pm
import pymaclab.modfiles.models as models

rbc = pm.newMOD(models.stable.rbc1_num,mesg=False,ncpus='auto')

# Try to update all of the wrapped objects and test if this has worked

# Do for paramdic, set_item
def test_paramdic_item():
    eta_key = 'eta'
    eta_old = 2.0
    eta_new = 5.0
    rbc.updaters.paramdic[eta_key] = eta_new
    
    # Did it work?
    assert rbc.updaters.paramdic.wrapobj[eta_key] == eta_new
    assert rbc.updaters.paramdic[eta_key] == eta_new
    assert rbc.paramdic[eta_key] == eta_new
    

# Do for paramdic, update
def test_paramdic_update():
    eta_key = 'eta'
    eta_old = 2.0
    eta_new = 5.0
    rho_key = 'rho'
    rho_old = 0.36
    rho_new = 0.35
    tmp_dic = {}
    tmp_dic[eta_key] = eta_new
    tmp_dic[rho_key] = rho_new
    rbc.updaters.paramdic.update(tmp_dic)
    
    # Did it work?
    assert rbc.updaters.paramdic.wrapobj[eta_key] == eta_new
    assert rbc.updaters.paramdic[eta_key] == eta_new
    assert rbc.paramdic[eta_key] == eta_new
    
    assert rbc.updaters.paramdic.wrapobj[rho_key] == rho_new
    assert rbc.updaters.paramdic[rho_key] == rho_new
    assert rbc.paramdic[rho_key] == rho_new
    
    

# Do for nlsubsdic, set_item
def test_nlsubsdic_item():
    U_key = '@U(t)'
    U_old = 'c(t)**(1-eta)/(1-eta)'
    U_new = 'c(t)**(1-eta*1.01)/(1-eta*1.01)'
    rbc.updaters.nlsubsdic[U_key] = U_new
    
    # Did it work?
    assert rbc.updaters.nlsubsdic.wrapobj[U_key] == U_new
    assert rbc.updaters.nlsubsdic[U_key] == U_new
    assert rbc.nlsubsdic[U_key] == U_new
    
    

# Do for nlsubsdic, update
def test_nlsubsdic_update():
    U_key = '@U(t)'
    U_old = 'c(t)**(1-eta)/(1-eta)'
    U_new = 'c(t)**(1-eta*1.01)/(1-eta*1.01)'
    F_key = '@F(t)'
    F_old = 'z(t)*k(t-1)**rho'
    F_new = 'z(t)*k(t-1)**rho*1.01'
    tmp_dic = {}
    tmp_dic[U_key] = U_new
    tmp_dic[F_key] = F_new
    rbc.updaters.nlsubsdic.update(tmp_dic)
    
    # Did it work?
    assert rbc.updaters.nlsubsdic.wrapobj[U_key] == U_new
    assert rbc.updaters.nlsubsdic[U_key] == U_new
    assert rbc.nlsubsdic[U_key] == U_new
    
    assert rbc.updaters.nlsubsdic.wrapobj[F_key] == F_new
    assert rbc.updaters.nlsubsdic[F_key] == F_new
    assert rbc.nlsubsdic[F_key] == F_new
    
    
    
# Do for vardic, set_item
def test_vardic_item():
    var_key = ['c(t)','consumption']
    indexo = rbc.vardic['con']['var'].index(var_key)
    var_old = 'bk'
    var_new = 'cf'
    rbc.updaters.vardic['con']['mod'][indexo][1] = var_new
    
    # Did it work?
    assert rbc.updaters.vardic.wrapobj['con']['mod'][indexo][1] == var_new
    assert rbc.updaters.vardic['con']['mod'][indexo][1] == var_new
    assert rbc.vardic['con']['mod'][indexo][1] == var_new
    
    