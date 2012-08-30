import os
from numpy import matlib as MAT
import numpy as np
import scikits.timeseries as ts
from scipy import io
import pymaclab as pm

#Specialize general options
#mac.mk_hessian = True
#mac.ncpus = 2
#NOTE: ?? use global flags or get rid of for a config file and automatization

#Add a datapath
datapath = '../data/'
# Add modpath
modpath = '../pymaclab/modfiles/'

# Get Christiano Data and prepare
input = open(os.path.join(datapath,'cee2005data.asc'),'r')
lines = input.read()
lines = lines.splitlines()
lines = lines[:-2]
i1=0
for x in lines:
   lines[i1] = x.split()[:]
   i2=0
   for y in lines[i1]:
      lines[i1][i2] = eval(y)
      i2 = i2 + 1
   i1 = i1 + 1


# Make sure all series are array(1,x)
linesm = MAT.matrix(lines)
chr_gdp = linesm[:,0].A.flatten(1)
chr_gdp = np.exp(chr_gdp/100).flatten(1)
chr_gdp = np.log(chr_gdp).flatten(1)
chr_pi = linesm[:,1].A.flatten(1)
chr_cons = linesm[:,2].A.flatten(1)
chr_cons = np.exp(chr_cons/100).flatten(1)
chr_cons = np.log(chr_cons).flatten(1)
chr_inv = linesm[:,3].A.flatten(1)
chr_inv = np.exp(chr_inv/100).flatten(1)

# make capital from investment
chr_cap = np.zeros((chr_inv.shape[0],))
chr_cap[0] = chr_inv[0]
for x in xrange(chr_inv.shape[0]-1):
    chr_cap[x+1] = chr_inv[x+1]+chr_cap[x]

chr_cap = np.log(chr_cap).flatten(1)
chr_inv = np.log(chr_inv).flatten(1)
chr_rwage = linesm[:,4].A.flatten(1)
chr_rwage = np.exp(chr_rwage/100).flatten(1)
chr_rwage = np.log(chr_rwage).flatten(1)
chr_prod = linesm[:,5].A.flatten(1)
chr_prod = np.exp(chr_prod/100).flatten(1)
chr_prod = np.log(chr_prod).flatten(1)
chr_ffr = linesm[:,6].A.flatten(1)
chr_dm = linesm[:,7].A.flatten(1)
chr_rprof = linesm[:,8].A.flatten(1)
chr_rprof = np.exp(chr_rprof/100).flatten(1)
chr_rprof = np.log(chr_rprof).flatten(1)


#chrt_err = ts.time_series\
         #(err,start_date=TS.Date(freq='Q', year=1964, quarter=2),freq='Q')
chrt_gdp = ts.time_series\
         (chr_gdp,start_date=ts.Date(freq='Q', year=1964, quarter=2),freq='Q')
chrt_pi = ts.time_series\
        (chr_pi,start_date=ts.Date(freq='Q', year=1964, quarter=2),freq='Q')
chrt_cons = ts.time_series\
        (chr_cons,start_date=ts.Date(freq='Q', year=1964, quarter=2),freq='Q')
chrt_inv = ts.time_series\
        (chr_inv,start_date=ts.Date(freq='Q', year=1964, quarter=2),freq='Q')
chrt_rwage = ts.time_series\
        (chr_rwage,start_date=ts.Date(freq='Q', year=1964, quarter=2),freq='Q')
chrt_prod = ts.time_series\
        (chr_prod,start_date=ts.Date(freq='Q', year=1964, quarter=2),freq='Q')
chrt_ffr = ts.time_series\
        (chr_ffr,start_date=ts.Date(freq='Q', year=1964, quarter=2),freq='Q')
chrt_dm = ts.time_series\
        (chr_dm,start_date=ts.Date(freq='Q', year=1964, quarter=2),freq='Q')
chrt_rprof = ts.time_series\
        (chr_rprof,start_date=ts.Date(freq='Q', year=1964, quarter=2),freq='Q')
chrt_cap = ts.time_series\
        (chr_cap,start_date=ts.Date(freq='Q', year=1964, quarter=2),freq='Q')

# Load Christiano's VAR data in matlab format
cdata = io.loadmat(os.path.join(datapath,'datrep.mat'))
cdata = cdata['datrep']
cdata = np.matrix(cdata)

# Some VAR Stuff
Var1 = pm.newVAR(4,cdata[:,:],'const')
Var1.ols()
Var1.ols_comp()
Var1.do_irf(3,15)


# Build database and populate
db1 = pm.newDB()
db1.nameimDatStr(os.path.join(datapath,'USEURSN.csv'))
#db1.tsim('chr_err','Steady-State error',chrt_err)
#db1.mkalias('chr_err','err')
db1.tsim('chr_gdp','GDP series from Christiano',chrt_gdp)
db1.mkhpf('chr_gdp','chr_gdp_f')
db1.mkalias('chr_gdp','output')
db1.tsim('chr_pi','Inflation Series from Christiano',chrt_pi)
db1.tsim('chr_cons','Consumption Series from Christiano',chrt_cons)
db1.mkhpf('chr_cons','chr_cons_f')
db1.mkalias('chr_cons','consumption')
db1.tsim('chr_inv','Investment Series from Christiano',chrt_inv)
db1.mkhpf('chr_inv','chr_inv_f')
db1.mkalias('chr_inv','capital')
db1.tsim('chr_rwage','Real Wage Series from Christiano',chrt_rwage)
db1.mkhpf('chr_rwage','chr_rwage_f')
db1.tsim('chr_prod','Productivity Series from Christiano',chrt_prod)
db1.tsim('chr_ffr','Real Fed Fund Rate Series from Christiano',chrt_ffr)
db1.mkalias('chr_ffr','interest')
db1.tsim('chr_dm','Money Supply Change Series from Christiano',chrt_dm)
db1.tsim('chr_rprof','Real Profit Series from Christiano',chrt_rprof)
db1.tsim('chr_cap','Capital Stock Series from Christiano',chrt_cap)
db1.mkhpf('chr_cap','chr_cap_f')
db1.mkmodif('Q')
varord = [['output',2],['consumption',1],['interest',3]]
varlag = 2
shock_pos = 2
#TODO This does not seem to work at the moment, db has no VAR instance error thrown...
#db1.mkvar(varlag,varord,shock_pos,'const')


# Some dsge stuff
#TODO Let's not make calls to the DSGE model with a dataset at the moment as it is broken
#rbc1 = pm.newMOD(os.path.join(datapath,'rbc1.txt'),db1)

rbc1 = pm.newMOD(os.path.join(modpath,'rbc1.txt'))
rbc1.ccv('forkleind')
sims = pm.newMOD(os.path.join(modpath,'sims.txt'))
sims.ccv('forkleind')
rbc2 = pm.newMOD(os.path.join(modpath,'rbc2.txt'))
rbc2.ccv('forkleind')
'''
mbc1 = pm.newMOD(os.path.join(modpath,'mbc1.txt'))
mbc1.ccv('forkleind')
model2 = pm.newMOD(os.path.join(modpath,'model2.txt'))
model2.ccv('forkleind')
model3 = pm.newMOD(os.path.join(modpath,'model3.txt'))
model3.ccv('forkleind')
#cmod = newMOD('max1.txt')
#cmod.ccv('forkleind')
cmod2 = pm.newMOD(os.path.join(modpath,'max2.txt'))
cmod2.ccv('forkleind')
nknc = pm.newMOD(os.path.join(modpath,'nk_nocapital.txt'))
nknc.ccv('forkleind')
'''
merz = pm.newMOD(os.path.join(modpath,'merz.txt'),initlev=1)
#merz.ccv('forkleind')
'''
cee = pm.newMOD(os.path.join(modpath,'cee.txt'))
cee.ccv('forkleind')
'''