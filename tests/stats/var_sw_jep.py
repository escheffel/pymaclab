import pandas
import numpy
import matplotlib
from matplotlib import pyplot as plt
from pymaclab import newVAR
import os

modconfig = {}
# This model's name or title, used for pickling results
modname = 'STOCKWATSON_jep_model'
modconfig['modname'] = modname
# Get the data from the directory
datafile = 'stwat_jpe_2001.csv'
modconfig['datafile'] = datafile
datapath = '../data/'
modconfig['datapath'] = datapath
data = pandas.read_csv(os.path.join(datapath,datafile),delimiter=';')
data = data[52:-1] # Select data range
data = data[:164] # Further restriction on data, because commodity index is longer
# Define frequency of data
freq = 'Q'
modconfig['freq'] = freq
# Should the data be standardized?
trans_data = False
modconfig['trans_data'] = trans_data
# Define the number of lags in the VAR
nlags = 4
modconfig['nlags'] = nlags
# How many periods for IRFs ?
irfp = 25
modconfig['irfp'] = irfp
# Calculate variance decomposition?
do_vdc = True
modconfig['do_vdc'] = do_vdc
# How many periods forecast for variance decomposition ?
vdcp = 20
modconfig['vdcp'] = vdcp
# Use the svname matrix for alternative shocks?
use_svnames = False
modconfig['use_svnames'] = use_svnames
# Transform log variables in IRFs into percentages ?
translog = False
modconfig['translog'] = translog
# Transform the percentage responses back into absolutes from the mean ?
transperc = False
modconfig['transperc'] = transperc
# Add impact matrix to Cholesky responses ?
cholimp = False
modconfig['cholimp'] = cholimp
# How many matrices for MA Representation
maxphi = irfp+2
modconfig['maxphi'] = maxphi
# Boostrap IRF confidence intervals ?
bootstrap = True
modconfig['bootstrap'] = bootstrap
####################### This is experimental stuff, need to clarify the difference between the two #####################
killian_bsinbs = True
modconfig['killian_bsinbs'] = killian_bsinbs
########################################################################################################################
########################################## Graphing options for IRFs and other graphs ##################################
graph_options = {}
graph_options['labels'] = {}
graph_options['labels']['title'] = True
graph_options['labels']['xlabel'] = False
graph_options['labels']['ylabel'] = True
graph_options['labels']['text_fs'] = 10
graph_options['labels']['ylabel_fs'] = 10
graph_options['labels']['xlabel_fs'] = 10
graph_options['labels']['title_fs'] = 12
graph_options['labels']['label_fs'] = 8
graph_options['labels']['xticks_fs'] = 8
graph_options['labels']['yticks_fs'] = 8
graph_options['labels']['legend_fs'] = 8
graph_options['labels']['use_tex'] = False
graph_options['lines'] = {}
graph_options['lines']['width'] = 2.0
graph_options['grid'] = True
graph_options['irfs'] = {}
graph_options['irfs']['outer_colour'] = '#D3B0B0'
graph_options['irfs']['inner_colour'] = '#EED6D6'
graph_options['irfs']['line_colour'] = 'black'
graph_options['irfs']['biased_line'] = True # This allows you to also show the biased IRF when Killian's b-i-b option is on
graph_options['irfs']['biased_line_type'] = 'g--'
graph_options['save'] = {}
graph_options['save']['format'] = {}
graph_options['save']['format']['eps'] = True
graph_options['save']['format']['pdf'] = True
modconfig['graph_options'] = graph_options
#########################################################################################################################
# Draw a horizontal zero line into each IRF Plot
hzline = True
modconfig['hzline'] = hzline
# How thick should the hline be?
hzline_width = 1.0
modconfig['hzline_width'] = hzline_width
# How many draws for bootstrap ?
bdraw = 1000
modconfig['bdraw'] = bdraw
# Use multicore processors to speed up bootstrap?
multicpu = True
modconfig['multicpu'] = multicpu
# Specify number of cpu cores or set to 'auto'
ncpus = 'auto'
modconfig['ncpus'] = ncpus
# Level of significance for bootstrap confidence intervals
signif = [0.66,0.95] # implies 1-alpha = 0.66, use 'mix' to graph a variety of CI
modconfig['signif'] = signif
# What kind of intercept/time trend?
const_term = 'ct' # Can be 'None', 'cc' or 'tt'
modconfig['const_term'] = const_term
# Pickle the model run results and settings
pickle_switch = False
# Define the names of the vars used
vnames = []
vnames.append('INF')
vnames.append('LHUR')
vnames.append('FYFF')
modconfig['vnames'] = vnames
# Define the names for the vnames used in plotting
pnames = {}
pnames['INF'] = 'Inflation'
pnames['LHUR'] = 'Unemployment'
pnames['FYFF'] = 'Federal Funds Rate'

# Create the data modification dictionary
tcode_dic = {}
clist = [1,1,1]
for i1,vname in enumerate(vnames):
    tcode_dic[vname] = clist[i1]
modconfig['tcode'] = tcode_dic

# Define a matrix for scaling the orthogonalized shock
svnames = {}
for name in vnames:
    svnames[name] = 1.0
modconfig['svnames'] = svnames

# Set the transperc_dic
transperc_dic = {}
for name in vnames:
    transperc_dic[name] = False
transperc_dic['employment'] = True
modconfig['transperc_dic'] = transperc_dic

# Get the individual series of interest
loglist = []
modconfig['loglist'] = loglist
tlist = []
for name in vnames:
    if name in loglist:
        tlist.append(data[name].values)
    else:
        tlist.append(data[name].values)
matd = numpy.array(tlist)
# Transpose and get the dims
matd = matd.T
cols = matd.shape[1]
rows = matd.shape[0]
modconfig['cols'] = cols
modconfig['rows'] = rows

swvar = newVAR(data=matd,vnames=vnames,pnames=pnames,svnames=svnames,irfs=True,boot=True,plot=True,conf=modconfig,mesg=True)
'''
tsvar = var_model.VAR(matd,['INF','LHUR','FYFF'])
tsvar_res = tsvar.fit(4,trend='ct')
'''

