import pymaclab as pm
import pymaclab.modfiles.models as models

# Define the ssidic of initial guesses or starting values
alpha = 0.3
theta = 0.85
gamma = 0.25
delta_k = 0.025
betta = 0.99
B     = 0.65
inf_bar = 1.012
ssidic = {}
ssidic['k_bar']    = 30.0
ssidic['n_bar']    = 0.33
ssidic['inf_bar']  = inf_bar
ssidic['r_bar']    = 1/betta
ssidic['y_bar']    = ssidic['k_bar']**alpha*ssidic['n_bar']**(1-alpha)
ssidic['c_bar']    = ssidic['y_bar'] - delta_k*ssidic['k_bar']
ssidic['m_bar']    = 5.0
# Define marginal cost
rk = alpha*(ssidic['n_bar']/ssidic['k_bar'])**(1-alpha)
rn = (1-alpha)*(ssidic['k_bar']/ssidic['n_bar'])**alpha
mc = (1/(1-alpha))**(1-alpha)*(1/alpha)**alpha*(rk)**alpha*(rn)
ssidic['aa_bar']   = inf_bar*mc*ssidic['y_bar']/(1-(1-gamma)*betta*ssidic['inf_bar']**((2-theta)/(1-theta)))
ssidic['bb_bar']   = ssidic['y_bar']/(1-(1-gamma)*betta*ssidic['inf_bar']**(1/(1-theta)))
ssidic['infr_bar'] = (1/theta)*(ssidic['aa_bar']/ssidic['bb_bar'])
ssidic['i_bar']    = ssidic['r_bar']*inf_bar
ssidic['lam_bar']  = 1.0/ssidic['c_bar']
ssidic['w_bar'] = (1/theta)*rn

# Instantiate a new DSGE model instance like so
foc_li = [0,1,2,3,4,5,6,7,8,9,10,11,12]
nkm1 = pm.newMOD(models.nkm,use_focs=foc_li,ssidic=ssidic,initlev=0) 