import pandas
import numpy
import matplotlib
from matplotlib import pyplot as plt
from pymaclab import newFAVAR
import os
import copy
import datetime as dt
from scikits import timeseries as ts
from scikits.timeseries.lib.interpolate import interp_masked1d,backward_fill
from pymaclab.dattrans import fred
import glob
fred_key = 'c5c6ec1dd5e1eb81b2727053690986db'
fred.key(fred_key)

modconfig = {}
# This model's name or title, used for pickling results
modname = 'BLOOM_favar_model'
modconfig['modname'] = modname
# Get the data from the directory, later we will add the St. Louis Fed
datafile = 'October_Updated_finaldata.csv'
modconfig['datafile'] = datafile
datapath = '../data/'
modconfig['datapath'] = datapath
data = pandas.read_csv(os.path.join(datapath,datafile),delimiter=';')
data = pandas.DataFrame(numpy.array([data['uncert'].values,data['spindex'].values,data['ffr'].values,
                                     data['employment'].values,data['ip'].values,data['copper'].values]).T,
                        columns=['UNCERT','SPINDEX','FFR','EMP','IP','COPPER'])
columns = list(data.columns)
# Name list for St. Louis Fed data to be appended
stfedli = ['GS2','GS10','PPIFGS','CPIAUCSL','CPILFESL','M2SL','CFNAI','IPDCONGD',
           'TCU','MCUMFN','HOUST','NAPM','USSLIND','IC4WSA','PSAVERT','MZM','TB3MS',
           'TB6MS','GS1','GS5','UMCSENT','LNS12032195','AAA','EMRATIO','OILPRICE',
           'DSPIC96','MORTG','CURRENCY','BAA','PCE','TOTALSL','GPDI','NAPMII',
           'NAPMNOI','NAPMEI','NAPMSDI','NAPMPRI','NAPMPI','CPIUFDNS','CUSR0000SAC',
           'CUSR0000SEHF','CUSR0000SAS4','PPIITM','GS7','GS3','PCEPI','PCEPILFE',
           'PCEDG','PCEND','PCES','NONREVSL','TWEXMMTH','MZMV','M2V','M1V']

for elem in stfedli:
    columns.append(elem)
# Define frequency of data
freq = 'M'
modconfig['freq'] = freq
###################################################################################################################
######################## Define here the lags for the data and the factors ########################################
### There is a trade-off here, you can chose not to use any lags on the data itself and let the factors capture ###
### all of the data's persistence, in that case your lag structure on the factors needs to be long enough #########
# Should the original data be filtered using their own lags?
xdata_filter = False
modconfig['xdata_filter'] = xdata_filter
# Define the number of lags in the VAR, only applicable if option above was chosen to be true
nlags = 6
modconfig['nlags'] = nlags
# Define the number of lags in the VAR of the factors, we always keep this at 1, becaus we are using a companion structure
flags = 6
modconfig['flags'] = flags
# Define the number of common static factors to be employed, but this can later be overriden by automatic selection
sfacs = 12
modconfig['sfacs'] = sfacs
###################################################################################################################
###################################################################################################################
# Define the convergence criterion in PERCENTAGES
conv_crit = 0.001
modconfig['conv_crit'] = conv_crit
# Should IRFs and confidence intervals be computed?
do_irfs = True
modconfig['do_irfs'] = do_irfs
# Should all directories be deleted before a model run?
predel = True
modconfig['predel'] = predel
# How many periods for IRFs ?
irfp = 48
modconfig['irfp'] = irfp
# Calculate variance decomposition?
do_vdc = True
modconfig['do_vdc'] = do_vdc
# How many periods forecast for variance decomposition ?
vdcp = 60
modconfig['vdcp'] = vdcp
# Graph only IR with respect to identified variable
irf_ident_graph = True
modconfig['irf_ident_graph'] = irf_ident_graph
# How many periods for the XDATA IRFs ?
xirfp = 40
modconfig['xirfp'] = xirfp
# Use the svname matrix for alternative shocks?
use_svnames = False
modconfig['use_svnames'] = use_svnames
# Rescale the IRFs ?
rescale = True
modconfig['rescale'] = rescale
# Transform log variables in IRFs into percentages (only useful with rescaling) ?
translog = True
modconfig['translog'] = translog
# Transform the percentage responses back into absolutes from the mean ?
transperc = False
modconfig['transperc'] = transperc
# Add impact matrix to Cholesky responses ?
cholimp = False
modconfig['cholimp'] = cholimp
# How many matrices for MA Representation, this has to be a bit bigger than the number of time periods in the IRs
maxphi = max(vdcp,irfp,xirfp)+2
modconfig['maxphi'] = maxphi
# Boostrap IRF confidence intervals ?
bootstrap = True
modconfig['bootstrap'] = bootstrap
# Draw a horizontal zero line into each IRF Plot
hzline = True
modconfig['hzline'] = hzline
# Should IRFs be annotated
irf_anot = False
modconfig['irf_anot'] = irf_anot
# How thick should the hline be?
hzline_width = 1.0
modconfig['hzline_width'] = hzline_width
# How many draws for bootstrap ?
bdraw = 1000
modconfig['bdraw'] = bdraw
# Use multicore processors to speed up bootstrap?
multicpu = True
modconfig['multicpu'] = multicpu
# If bootstraps are split into chunks, what is each chunk's size (matters for memory efficiency)
bdraw_chunk = 200
modconfig['bdraw_chunk'] = bdraw_chunk
# Specify number of cpu cores or set to 'auto'
ncpus = 'auto'
modconfig['ncpus'] = ncpus
# Level of significance for bootstrap confidence intervals
signif = 'mix' # implies 1-alpha = 0.66, use 'mix' to graph a variety of CI
modconfig['signif'] = signif
# Boostrap with replacement ?
breplace = False
modconfig['breplace'] = breplace
# What kind of intercept/time trend?
const_term = None # Can be 'None', 'cc' or 'tt'
modconfig['const_term'] = const_term
# Pickle the model run results and settings
pickle_switch = False
modconfig['pickle_switch'] = pickle_switch
########################### use ONE of the identification schemes, not many ###################
# Employ BBE setup of the FAVAR?
use_bbe_ident = True
modconfig['use_bbe_ident'] = use_bbe_ident
# Employ simple Cholesky factorisation of the factor innovation matrix to diagonalize shocks
use_chol_ident = False
modconfig['use_chol_ident'] = use_chol_ident
################################################################################################
####################### This is experimental stuff, need to clarify the difference between the two #####################
# Parametric bootstrap vs. non-parametric bootstrap?
param_bootstrap = False
modconfig['param_bootstrap'] = param_bootstrap
killian_bsinbs = True
modconfig['killian_bsinbs'] = killian_bsinbs
########################################################################################################################
########################################## Graphing options for IRFs and other graphs ##################################
graph_options = {}
graph_options['labels'] = {}
graph_options['labels']['title'] = True
graph_options['labels']['title_top_only'] = True
graph_options['labels']['xlabel'] = False
graph_options['labels']['ylabel'] = True
graph_options['labels']['xlabel_fs'] = 10
graph_options['labels']['ylabel_fs'] = 10
graph_options['labels']['text_fs'] = 10
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
graph_options['irfs']['outer_colour'] = '#a4a4a4'
graph_options['irfs']['inner_colour'] = '#c6c6c6'
graph_options['irfs']['line_colour'] = 'black'
graph_options['irfs']['biased_line'] = True # This allows you to also show the biased IRF when Killian's b-i-b option is on
graph_options['irfs']['biased_line_type'] = 'k--'
modconfig['graph_options'] = graph_options
#########################################################################################################################
# Should recessions be shown in the factor graphs?
plot_rec = True
modconfig['plot_rec'] = plot_rec
# Colour and transparency of recession shades
rec_colour = '#a4a4a4'
rec_transp = 0.2
modconfig['rec_colour'] = rec_colour
modconfig['rec_transp'] = rec_transp
##########################################################################################################################
# Options for BBE identification setup
separate_vars = []
separate_vars.append('UNCERT')
modconfig['separate_vars'] = separate_vars
slow_vars = []
slow_vars.append('EMP')
slow_vars.append('IP')
slow_vars.append('IPDCONGD')
slow_vars.append('PPIFGS')
slow_vars.append('TCU')
slow_vars.append('MCUMFN')
slow_vars.append('CPIAUCSL')
slow_vars.append('CPILFESL')
slow_vars.append('HOUST')
slow_vars.append('IC4WSA')
slow_vars.append('PSAVERT')
slow_vars.append('LNS12032195')
slow_vars.append('EMRATIO')
slow_vars.append('DSPIC96')
slow_vars.append('PCE')
slow_vars.append('GPDI')
slow_vars.append('CPIUFDNS')
slow_vars.append('CUSR0000SAC')
slow_vars.append('CUSR0000SEHF')
slow_vars.append('CUSR0000SAS4')
slow_vars.append('PPIITM')
slow_vars.append('PCEPI')
slow_vars.append('PCEPILFE')
slow_vars.append('PCEDG')
slow_vars.append('PCEND')
slow_vars.append('PCES')
slow_vars.append('MZMV')
slow_vars.append('M2V')
slow_vars.append('M1V')
modconfig['slow_vars'] = slow_vars
fast_vars = []
fast_vars.append('COPPER')
fast_vars.append('OILPRICE')
fast_vars.append('CFNAI')
fast_vars.append('NAPM')
fast_vars.append('NAPMII')
fast_vars.append('NAPMNOI')
fast_vars.append('NAPMEI')
fast_vars.append('NAPMSDI')
fast_vars.append('NAPMPRI')
fast_vars.append('NAPMPI')
fast_vars.append('UMCSENT')
fast_vars.append('NONREVSL')
fast_vars.append('SPINDEX')
fast_vars.append('FFR')
fast_vars.append('GS10')
fast_vars.append('GS2')
fast_vars.append('M2SL')
fast_vars.append('MZM')
fast_vars.append('USSLIND')
fast_vars.append('TB3MS')
fast_vars.append('TB6MS')
fast_vars.append('GS1')
fast_vars.append('GS5')
fast_vars.append('AAA')
fast_vars.append('MORTG')
fast_vars.append('CURRENCY')
fast_vars.append('BAA')
fast_vars.append('TOTALSL')
fast_vars.append('GS7')
fast_vars.append('GS3')
fast_vars.append('TWEXMMTH')

modconfig['fast_vars'] = fast_vars
# You can also fix the number of factors for the slow variables here, or let them be estimated using Ng-Bai criterion
slow_facs = sfacs-2
modconfig['slow_facs'] = slow_facs
# Define the names of the vars used
vnames = []
for sname in columns:
    vnames.append(sname)
modconfig['vnames'] = vnames
# Define the alternative shocks
svnames = {}
for sname in columns:
    svnames[sname] = 1.0
svnames['UNCERT'] = 112.0
modconfig['svnames'] = vnames
# Set the transperc_dic
transperc_dic = {}
for name in vnames:
    transperc_dic[name] = False
modconfig['transperc_dic'] = transperc_dic
# Define the names for the vnames used in plotting, at the moment copy the short name
pnames = {}
for sname in columns:
    pnames[sname] = sname
#Define the variables which are cointegrated, here the order also matters !
conames = []
conames.append('INF')
conames.append('LHUR')
conames.append('FYFF')
modconfig['conames'] = conames
# Define the matrix betta for the long-run cointegrating relationship
betta_coint = numpy.array([[1,2,3],[1,2,3],[1,2,3]])
modconfig['betta_coint'] = betta_coint


# Function to transform the raw data according to some classification
# 1 = levels
# 2 = first seasonal difference
# 3 = second seasonal difference
# 4 = log level
# 5 = log first seasonal difference
# 6 = log second seasonal difference
# 7 = detrend log using hp filter monthly data
# 8 = detrend log using hp filter quarterly data
# 16 = log second seasonal difference
# 17 = (1-L)(1-L^12)
# 18 = log of (1-L)*annualizing factor (i.e. x4 for quarterly and x12 for monthly
# 19 = filter the level using simple time-trend
# 20 = filter the log using simple time-trend
# Create the data modification dictionary
tcode_dic = {}
clist = [1,4,1,7,7,4]
stlouislist = [1,1,7,7,7,7,1,7,1,1,4,1,1,4,1,7,1,1,1,1,1,7,1,4,4,
               7,1,7,1,7,7,7,1,1,1,1,1,1,7,7,7,7,7,1,1,7,7,7,7,7,7,4,4,4,4]
clist = clist + stlouislist
for i1,vname in enumerate(vnames):
    tcode_dic[vname] = clist[i1]
modconfig['tcode'] = tcode_dic

# Define a matrix for scaling the orthogonalized shock
svnames = {}
for name in vnames:
    svnames[name] = 1.0
svnames['UNCERT'] = 112.0
modconfig['svnames'] = svnames
# Get the individual series of interest, also take natural logs where wanted
loglist = []
modconfig['loglist'] = loglist
orig_li = vnames[:-len(stfedli)]
tlist = []
for name in orig_li:
    if name in loglist:
        tlist.append(numpy.log(data[name].values))
    else:
        tlist.append(data[name].values)
matd = numpy.array(tlist)
matd = matd[:,:-2]
# Transpose and get the , but first add St. Lous Fed data, also save in temporary file
matd = matd.T
# Also create the table for the data
tmpstr = ''
tmpstr = tmpstr + r"\begin{sidewaystable}[H]\centering"+"\n"
tmpstr = tmpstr + r"\caption{Summary of all data used in FAVAR}"+"\n"
tmpstr = tmpstr + r"\tiny"+"\n"
tmpstr = tmpstr + r"\begin{tabular}{@{}lllllll@{}}\toprule"+"\n"
tmpstr = tmpstr + r"No. & Description & Code & Source & Original Units & Variable Type & TCode\\\midrule"+"\n"
newstr = str(1)+r" & "+r"Bloom et al.'s policy uncertainty index"+r" & "+r"UNCERT"+r" & "+r"Bloom et al."+r" & "+r"Index"+r" & "+r"Identified"+r" & "+str(4)+r"\\"+"\n"
tmpstr = tmpstr + newstr
newstr = str(2)+r" & "+r"Standard and Poor's 500 US Stock Market Index"+r" & "+r"SPINDEX"+r" & "+r"Bloom et al."+r" & "+r"Index"+r" & "+r"Fast"+r" & "+str(4)+r"\\"+"\n"
tmpstr = tmpstr + newstr
newstr = str(3)+r" & "+r"US Federal Funds Rate"+r" & "+r"FFR"+r" & "+r"Bloom et al."+r" & "+r"Percent"+r" & "+r"Fast"+r" & "+str(1)+r"\\"+"\n"
tmpstr = tmpstr + newstr
newstr = str(4)+r" & "+r"US Total Civilian Employment"+r" & "+r"EMP"+r" & "+r"Bloom et al."+r" & "+r"Total no. of persons in 1000s"+r" & "+r"Slow"+r" & "+str(7)+r"\\"+"\n"
tmpstr = tmpstr + newstr
newstr = str(5)+r" & "+r"US Industrial Production Index"+r" & "+r"IP"+r" & "+r"Bloom et al."+r" & "+r"Index"+r" & "+r"Slow"+r" & "+str(7)+r"\\"+"\n"
tmpstr = tmpstr + newstr
newstr = str(6)+r" & "+r"Price of Copper in US dollars"+r" & "+r"COPPER"+r" & "+r"Bloom et al."+r" & "+r"US dollars per tonne"+r" & "+r"Fast"+r" & "+str(4)+r"\\"+"\n"
tmpstr = tmpstr + newstr
for i1,namo in enumerate(stfedli):
    print "Now acquiring data for: "+namo
    resu = fred.observations(namo)
    resuf = fred.series(series_id=namo)
    sera = resuf['seriess']['series']
    if namo in slow_vars: stringo = 'Slow'
    elif namo in fast_vars: stringo = 'Fast'
    elif namo in separate_vars: stringo = 'Identified'
    if namo == 'LNS12032195': stitle = sera['title'].split(',')[0]
    else: stitle = sera['title'].replace(r'&','and')
    newstr = str(i1+7)+r" & "+stitle+r" & "+namo+r" & "+r"St. Louis Fed"+r" & "+ sera['units']+r" & "+stringo+r" & "+str(stlouislist[i1])+r"\\"+"\n"
    tmpstr = tmpstr + newstr
    freqf = resuf['seriess']['series']['frequency_short']
    print "Done..."
    valo = [x['value'] for x in resu['observations']['observation']]
    valo2 = numpy.zeros(len(valo))
    for i1,elem in enumerate(valo):
        if elem == '.': valo2[i1] = float(valo2[i1-1])
        else: valo2[i1] = float(elem)
    valo = numpy.array(valo2)
    dato = [x['date'] for x in resu['observations']['observation']]
    if freqf == 'D':
        tso = ts.time_series(valo,dates=dato,freq='D')
        tso2 = ts.convert(tso,'M',func=numpy.ma.mean)
    elif freqf == 'W':
        tso = ts.time_series(valo,dates=dato,freq='D')
        tso2 = ts.convert(tso,'M',func=numpy.ma.mean)
    elif freqf == 'M':
        tso = ts.time_series(valo,dates=dato,freq='M')
        tso2 = copy.deepcopy(tso)
    elif freqf == 'Q':
        tso = ts.time_series(valo,dates=dato,freq='Q')
        tso2 = ts.convert(tso,'M')
        tso2 = interp_masked1d(tso2)
    valor = tso2['1985-01-01':'2011-09-01']
    dator = valor.dates
    valor = valor.data
    valor = numpy.reshape(valor,(len(valor),1))
    matd = numpy.hstack((matd,valor))

tmpstr = tmpstr + r"\bottomrule"+"\n"
tmpstr = tmpstr + r"\multicolumn{7}{l}{{\bf Notes:} Transformation codes are: 1=levels, 2=first seasonal difference, 3=second seasonal difference, 4=log level, 5=log first seasonal difference, 6=log second seasonal difference,}\\"+"\n"
tmpstr = tmpstr + r"\multicolumn{7}{l}{7=log hp-filtered monthly data. Note that investment data were intrapolated in order to obtain monthly from quarterly data. The choice of fast/slow variables in this table reflects the benchmark ordering.}"+"\n"
tmpstr = tmpstr + r"\end{tabular}"+"\n"
tmpstr = tmpstr + r"\end{sidewaystable}"+"\n"

'''
if 'table1.tex' in os.listdir('../tex_papers/bloom_favar/'):
    os.remove('../tex_papers/bloom_favar/table1.tex')
f1 = open('../tex_papers/bloom_favar/table1.tex','w')
f1.write(tmpstr)
f1.close()
'''

cols = matd.shape[1]
rows = matd.shape[0]
modconfig['cols'] = cols
modconfig['rows'] = rows
# In total there are three US recessions in the sample considered here
rdates = []
rdates.append([ts.Date(freq='M',year=1990,month=7,day=1),ts.Date(freq='M',year=1991,month=3,day=1)])
rdates.append([ts.Date(freq='M',year=2001,month=3,day=1),ts.Date(freq='M',year=2001,month=11,day=1)])
rdates.append([ts.Date(freq='M',year=2007,month=12,day=1),ts.Date(freq='M',year=2009,month=7,day=1)])

# Before running the VAR, plot the uncertainty data for the paper
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.dates as mdates
majorLocator   = MultipleLocator(1)
majorFormatter = FormatStrFormatter('%d')
minorLocator   = MultipleLocator(1)
years    = mdates.YearLocator()   # every year
months   = mdates.MonthLocator()  # every month
yearsFmt = mdates.DateFormatter('%Y') 
glopt = modconfig['graph_options']['labels']
gridopt = modconfig['graph_options']['grid']        
modconfig['modname']
# Check if directory exists
if modname not in os.listdir('../graphs/'):
    os.mkdir('../graphs/'+modname)
uncert = matd[:,0]
fig1 = plt.figure()
ax1 = fig1.add_subplot(1,1,1)
if gridopt: ax1.grid()
for label in ax1.get_xticklabels():
    label.set_rotation(55)       
if glopt['xlabel']: ax1.set_xlabel('Time',fontsize=glopt['label_fs'])
if glopt['ylabel']: ax1.set_ylabel('Index 100=(1985-2000)',fontsize=glopt['label_fs'])
datli = [x.datetime for x in dator]
ax1.plot(datli[-len(uncert):],uncert,'k-',linewidth=2.0)
dato1 = dt.datetime(1985,12,31)
if glopt['title']: ax1.set_title('Uncertainty Index (Bloom et al.)',fontsize=glopt['title_fs'])
fig1.savefig('../graphs/'+modname+'/'+'UNCERT.eps',bbox_inches='tight')

favar = newFAVAR(data=matd,dates=dator,rdates=rdates,freq=freq,vnames=vnames,pnames=pnames,svnames=svnames,irfs=True,rescale=rescale,boot=True,plot=True,sfacs=sfacs,init=None,conf=modconfig,mesg=True)


# Also create the table for the impulse responses
tmpstr = ''
tmpstr = tmpstr + r"\begin{tabular}{@{}lll}"+"\n"
substr = ''
rowo = 7
colo = 3
rowo_count = 0
switcho = False
globli = glob.glob('../graphs/'+modname+'/shock_irfs/irf_12_*.eps')
globli.sort()
for i1,filo in enumerate(globli):
    if i1 != 0 and float(rowo_count+1)%float(rowo) == 0.0:
        tmpstr = tmpstr + r"\end{tabular}"+"\n"+r"\newpage"+"\n"+r"\begin{tabular}{@{}lll}"+"\n"
        rowo_count += 1
        switcho = True
    if i1 != 0 and float(i1+1)%float(colo) > 0.0: substr = substr + r" & \includegraphics[scale=0.27]{"+"../"+filo[:-4]+"}"
    elif i1 == 0: substr = substr + r" & \includegraphics[scale=0.27]{"+"../"+filo[:-4]+"}"
    else:
        substr = substr + r" & \includegraphics[scale=0.27]{"+"../"+filo[:-4]+"}"
        if i1 != 0: tmpstr = tmpstr + copy.deepcopy(substr[2:] + r"\\" + "\n")
        if not switcho: rowo_count += 1
        substr = ''
    switcho = False
tmpstr = tmpstr + r"\end{tabular}"+"\n"
if 'table2.tex' in os.listdir('../tex_papers/bloom_favar/'):
    os.remove('../tex_papers/bloom_favar/table2.tex')
f1 = open('../tex_papers/bloom_favar/table2.tex','w')
f1.write(tmpstr)
f1.close()

# Also create the table for the impulse responses, copper and oil
tmpstr = ''
tmpstr = tmpstr + r"\begin{tabular}{@{}lll}"+"\n"
substr = ''
rowo = 7
colo = 3
rowo_count = 0
switcho = False
globli = glob.glob('../graphs/'+modname+'/oilcopper/irf_12_*.eps')
globli.sort()
for i1,filo in enumerate(globli):
    if i1 != 0 and float(rowo_count+1)%float(rowo) == 0.0:
        tmpstr = tmpstr + r"\end{tabular}"+"\n"+r"\newpage"+"\n"+r"\begin{tabular}{@{}lll}"+"\n"
        rowo_count += 1
        switcho = True
    if i1 != 0 and float(i1+1)%float(colo) > 0.0: substr = substr + r" & \includegraphics[scale=0.27]{"+"../"+filo[:-4]+"}"
    elif i1 == 0: substr = substr + r" & \includegraphics[scale=0.27]{"+"../"+filo[:-4]+"}"
    else:
        substr = substr + r" & \includegraphics[scale=0.27]{"+"../"+filo[:-4]+"}"
        if i1 != 0: tmpstr = tmpstr + copy.deepcopy(substr[2:] + r"\\" + "\n")
        if not switcho: rowo_count += 1
        substr = ''
    switcho = False
tmpstr = tmpstr + r"\end{tabular}"+"\n"
if 'table3.tex' in os.listdir('../tex_papers/bloom_favar/'):
    os.remove('../tex_papers/bloom_favar/table3.tex')
f1 = open('../tex_papers/bloom_favar/table3.tex','w')
f1.write(tmpstr)
f1.close()

# Also create the table for the impulse responses, copper oil survey
tmpstr = ''
tmpstr = tmpstr + r"\begin{tabular}{@{}lll}"+"\n"
substr = ''
rowo = 7
colo = 3
rowo_count = 0
switcho = False
globli = glob.glob('../graphs/'+modname+'/oilcoppersurvey/irf_12_*.eps')
globli.sort()
for i1,filo in enumerate(globli):
    if i1 != 0 and float(rowo_count+1)%float(rowo) == 0.0:
        tmpstr = tmpstr + r"\end{tabular}"+"\n"+r"\newpage"+"\n"+r"\begin{tabular}{@{}lll}"+"\n"
        rowo_count += 1
        switcho = True
    if i1 != 0 and float(i1+1)%float(colo) > 0.0: substr = substr + r" & \includegraphics[scale=0.27]{"+"../"+filo[:-4]+"}"
    elif i1 == 0: substr = substr + r" & \includegraphics[scale=0.27]{"+"../"+filo[:-4]+"}"
    else:
        substr = substr + r" & \includegraphics[scale=0.27]{"+"../"+filo[:-4]+"}"
        if i1 != 0: tmpstr = tmpstr + copy.deepcopy(substr[2:] + r"\\" + "\n")
        if not switcho: rowo_count += 1
        substr = ''
    switcho = False
tmpstr = tmpstr + r"\end{tabular}"+"\n"
if 'table4.tex' in os.listdir('../tex_papers/bloom_favar/'):
    os.remove('../tex_papers/bloom_favar/table4.tex')
f1 = open('../tex_papers/bloom_favar/table4.tex','w')
f1.write(tmpstr)
f1.close()


# Create commodity markets graph
fig2 = plt.figure()
bnames = ['COPPER','OILPRICE']
bcols = ['#000000','#B5B3B4']
pron = ['min','min']
btype = ['-','-']
indx = []
for namos in bnames:
    indx.append(favar.vnames.index(namos))
ax2 = fig2.add_subplot(1,1,1)
plt.grid()
for i1,elem in enumerate(indx):
    if tcode_dic[bnames[i1]] in [4,5,6,7,8]:
        percen = 100.0
    else:
        percen = 1.0
    ax2.plot(favar.phisc[:49,elem,-1].flatten()*percen,label=bnames[i1],color=bcols[i1],linestyle=btype[i1])
    ax2.set_ylabel('Percentage Deviation',fontsize=glopt['label_fs'])
    '''
    if pron[i1] == 'max':
        val = max(favar.phisc[:21,elem,-1].flatten())
    elif pron[i1] == 'min':
        val = min(favar.phisc[:21,elem,-1].flatten())        
    pos = list(favar.phisc[:21,elem,-1].flatten()).index(val)
    if tcode_dic[bnames[i1]] in [5,7]:
        anval = numpy.round(favar.trans_dic[bnames[i1]]['std']*val*100.0,3)
    else:
        anval = numpy.round(favar.trans_dic[bnames[i1]]['std']*val,3)        
    ax2.annotate(str(anval), xy=(pos, val), xytext=(-5,5), 
                textcoords='offset points', ha='center', va='bottom', color='k')
    '''
if glopt['title']: ax2.set_title(r'Commodity market effects from shock to political uncertainty',fontsize=glopt['title_fs'])
ax2.hlines(numpy.array([0,]*48),xmin=0,xmax=48,linestyles='solid',linewidth=hzline_width,colors='k')
plt.legend(loc="lower right")
fig2.savefig('../graphs/'+modname+'/'+'COMM.eps',bbox_inches='tight')


# Create flight to safety graph
fig2 = plt.figure()
bnames = ['TB3MS','TB6MS','GS1','GS2','GS3','GS5','GS7','GS10']
bcols = ['#000000','#302F30','#6B6B6B','#B5B3B4','#000000','#302F30','#6B6B6B','#B5B3B4']
pron = ['min','min','min','min','min','min','min','min']
btype = ['-','-','-','-','--','--','--','--']
indx = []
for namos in bnames:
    indx.append(favar.vnames.index(namos))
ax2 = fig2.add_subplot(1,1,1)
plt.grid()
for i1,elem in enumerate(indx):
    if tcode_dic[bnames[i1]] in [4,5,6,7,8]:
        percen = 100.0
    else:
        percen = 1.0
    ax2.plot(favar.phisc[:49,elem,-1].flatten()*percen,label=bnames[i1],color=bcols[i1],linestyle=btype[i1])
    ax2.set_ylabel('Percentage Deviation',fontsize=glopt['label_fs'])
    '''
    if pron[i1] == 'max':
        val = max(favar.phisc[:21,elem,-1].flatten())
    elif pron[i1] == 'min':
        val = min(favar.phisc[:21,elem,-1].flatten())        
    pos = list(favar.phisc[:21,elem,-1].flatten()).index(val)
    if tcode_dic[bnames[i1]] in [5,7]:
        anval = numpy.round(favar.trans_dic[bnames[i1]]['std']*val*100.0,3)
    else:
        anval = numpy.round(favar.trans_dic[bnames[i1]]['std']*val,3)        
    ax2.annotate(str(anval), xy=(pos, val), xytext=(-5,5), 
                textcoords='offset points', ha='center', va='bottom', color='k')
    '''
if glopt['title']: ax2.set_title(r'Bond market effects from shock to political uncertainty',fontsize=glopt['title_fs'])
ax2.hlines(numpy.array([0,]*48),xmin=0,xmax=48,linestyles='solid',linewidth=hzline_width,colors='k')
plt.legend(loc="lower right")
fig2.savefig('../graphs/'+modname+'/'+'SAFETY.eps',bbox_inches='tight')

# Create external finance premium graph
fig2 = plt.figure()
bnames = ['AAA','BAA','SPINDEX']
bcols = ['#000000','#302F30','#6B6B6B']
pron = ['max','min','min']
btype = ['--','--','-']
indx = []
for namos in bnames:
    indx.append(favar.vnames.index(namos))
ax2 = fig2.add_subplot(1,1,1)
ax21 = plt.twinx()
plt.grid()
linli = []
for i1,elem in enumerate(indx):
    if tcode_dic[bnames[i1]] in [4,5,6,7,8]:
        percen = 100.0
    else:
        percen = 1.0
    if bnames[i1] == 'SPINDEX':
        linli.append(ax2.plot(favar.phisc[:49,elem,-1].flatten()*percen,label=bnames[i1],color=bcols[i1],linestyle=btype[i1]))
    elif bnames[i1] == 'AAA' or bnames[i1] == 'BAA':
        linli.append(ax21.plot(favar.phisc[:49,elem,-1].flatten()*percen,label=bnames[i1],color=bcols[i1],linestyle=btype[i1]))
    ax2.set_ylabel('Percentage Deviation',fontsize=glopt['label_fs'])
    ax21.set_ylabel('Percentage Deviation',fontsize=glopt['label_fs'])
    '''
    if pron[i1] == 'max':
        val = max(favar.phisc[:21,elem,-1].flatten())
    elif pron[i1] == 'min':
        val = min(favar.phisc[:21,elem,-1].flatten())        
    pos = list(favar.phisc[:21,elem,-1].flatten()).index(val)
    if tcode_dic[bnames[i1]] in [5,7]:
        anval = numpy.round(favar.trans_dic[bnames[i1]]['std']*val*100.0,3)
    else:
        anval = numpy.round(favar.trans_dic[bnames[i1]]['std']*val,3)        
    ax2.annotate(str(anval), xy=(pos, val), xytext=(-5,5), 
                textcoords='offset points', ha='center', va='bottom', color='k')
    '''
if glopt['title']: ax2.set_title(r'External finance effect from shock to political uncertainty',fontsize=glopt['title_fs'])
ax2.hlines(numpy.array([0,]*48),xmin=0,xmax=48,linestyles='solid',linewidth=hzline_width,colors='k')
# added these three lines
labs = [l[0].get_label() for l in linli]
labs[0] = labs[0]+"      (right -->)"
labs[1] = labs[1]+"      (right -->)"
labs[2] = labs[2]+" (<-- left)"
ax2.legend([x[0] for x in linli], labs, loc="lower right")

fig2.savefig('../graphs/'+modname+'/'+'EXTFIN.eps',bbox_inches='tight')

# Create consumption-saving graph
fig2 = plt.figure()
bnames = ['PCES','PCEND','PCEDG','PCEPI','PSAVERT']
pron = ['min','min','min','min','max']
bcols = ['#000000','#302F30','#6B6B6B','#000000','#6B6B6B']
btype = ['-','-','-','--','--']
indx = []
for namos in bnames:
    indx.append(favar.vnames.index(namos))
ax2 = fig2.add_subplot(1,1,1)
plt.grid()
for i1,elem in enumerate(indx):
    if tcode_dic[bnames[i1]] in [4,5,6,7,8]:
        percen = 100.0
    else:
        percen = 1.0    
    ax2.plot(favar.phisc[:49,elem,-1].flatten()*percen,label=bnames[i1],color=bcols[i1],linestyle=btype[i1])
    ax2.set_ylabel('Percentage Deviation',fontsize=glopt['label_fs'])
    '''
    if pron[i1] == 'max':
        val = max(favar.phisc[:21,elem,-1].flatten())
    elif pron[i1] == 'min':
        val = min(favar.phisc[:21,elem,-1].flatten())        
    pos = list(favar.phisc[:21,elem,-1].flatten()).index(val)
    if tcode_dic[bnames[i1]] in [4,5,6,7,8]:
        anval = numpy.round(favar.trans_dic[bnames[i1]]['std']*val*100.0,3)
    else:
        anval = numpy.round(favar.trans_dic[bnames[i1]]['std']*val,3)        
    ax2.annotate(str(anval), xy=(pos, val), xytext=(-5,5), 
                textcoords='offset points', ha='center', va='bottom', color='k')
    '''
if glopt['title']: ax2.set_title(r'Consumption & saving effect from shock to political uncertainty',fontsize=glopt['title_fs'])
ax2.hlines(numpy.array([0,]*48),xmin=0,xmax=48,linestyles='solid',linewidth=hzline_width,colors='k')
plt.legend(loc="lower right")
fig2.savefig('../graphs/'+modname+'/'+'CONSSAV.eps',bbox_inches='tight')


fevdli = []
fevdli.append('UNCERT')
fevdli.append('IP')
fevdli.append('EMP')
fevdli.append('SPINDEX')
fevdli.append('PCES')
fevdli.append('PCEDG')
fevdli.append('PCEND')
fevdli.append('PCES')
fevdli.append('TB3MS')
fevdli.append('TB6MS')
fevdli.append('GS1')
fevdli.append('GS2')
fevdli.append('GS5')

fevdli = copy.deepcopy(vnames)

# Also create the table for the variance decomposition
tmpstr = ''
tmpstr = tmpstr + r"\begin{table}\centering"+"\n"
tmpstr = tmpstr + r"\caption{Contribution of Policy Uncertainty Shock to Variance of Forecasts\\Benchmark ordering}"+"\n"
tmpstr = tmpstr + r"\footnotesize"+"\n"
tmpstr = tmpstr + r"\begin{tabular*}{0.8\textwidth}{@{\extracolsep{\fill}}@{}llccccccc@{}}\hline\hline"+"\n"
tmpstr = tmpstr + r" & & \multicolumn{6}{c}{Variance Decomposition at h}"+r"\\"+"\n"
tmpstr = tmpstr + r"\cline{3-8}"+"\n"
tmpstr = tmpstr + r"No. & Variable & h=6 & h=12 & h=18 & h=24 & h=36 & h=60 & $R^2$"+r"\\"+"\n"
tmpstr = tmpstr + r"\hline"+"\n"
for i1,namo in enumerate(fevdli):
    tmpstr = tmpstr + str(i1) + r" & "
    tmpstr = tmpstr + namo
    val6 = numpy.round(favar.decomp[vnames.index(namo),5,-1]/100.0,3)
    val6s = str(val6)
    if len(val6s) == 4:
        val6s = val6s + '0'
    elif len(val6s) == 3:
        val6s = val6s + '00'    
    tmpstr = tmpstr + r"&"+ val6s
    val12 = numpy.round(favar.decomp[vnames.index(namo),11,-1]/100.0,3)
    val12s = str(val12)
    if len(val12s) == 4:
        val12s = val12s + '0'
    elif len(val12s) == 3:
        val12s = val12s + '00'    
    tmpstr = tmpstr + r"&"+ val12s    
    val18 = numpy.round(favar.decomp[vnames.index(namo),17,-1]/100.0,3)
    val18s = str(val18)
    if len(val18s) == 4:
        val18s = val18s + '0'
    elif len(val18s) == 3:
        val18s = val18s + '00'    
    tmpstr = tmpstr + r"&"+ val18s
    val24 = numpy.round(favar.decomp[vnames.index(namo),23,-1]/100.0,3)
    val24s = str(val24)
    if len(val24s) == 4:
        val24s = val24s + '0'
    elif len(val24s) == 3:
        val24s = val24s + '00'    
    tmpstr = tmpstr + r"&"+ val24s
    val36 = numpy.round(favar.decomp[vnames.index(namo),35,-1]/100.0,3)
    val36s = str(val36)
    if len(val36s) == 4:
        val36s = val36s + '0'
    elif len(val36s) == 3:
        val36s = val36s + '00'    
    tmpstr = tmpstr + r"&"+ val36s
    val60 = numpy.round(favar.decomp[vnames.index(namo),59,-1]/100.0,3)
    val60s = str(val60)
    if len(val60s) == 4:
        val60s = val60s + '0'
    elif len(val60s) == 3:
        val60s = val60s + '00'    
    tmpstr = tmpstr + r"&"+ val60s
    yy = favar.tdata_adj[:,vnames.index(namo)]
    xx = favar.fdata_adj
    betta = numpy.dot(numpy.linalg.inv(numpy.dot(xx.T,xx)),numpy.dot(xx.T,yy))
    yyfit = numpy.dot(xx,betta.T)
    RSS = numpy.sum((yy - yyfit)**2)
    TSS = numpy.sum((yy - numpy.mean(yy))**2)
    ESS = TSS - RSS
    R2 = numpy.round(ESS/TSS,3)
    sR2 = str(R2)
    if len(sR2) == 4:
        sR2 = sR2 + '0'
    elif len(sR2) == 3:
        sR2 = sR2 + '00'
    if namo != "UNCERT":
        tmpstr = tmpstr + r"&"+ sR2+r"\\"+"\n"
    else:
        tmpstr = tmpstr + r"&"+ "1.000"+r"\\"+"\n"
tmpstr = tmpstr + r"\hline\hline"+"\n"
tmpstr = tmpstr + r"\multicolumn{9}{l}{{\bf Notes:} Right-hand side column shows $R^2$ from regression of variables on factors.}\\"+"\n"
tmpstr = tmpstr + r"\end{tabular*}"+"\n"
tmpstr = tmpstr + r"\end{table}"+"\n"
if 'table5.tex' in os.listdir('../tex_papers/bloom_favar/'):
    os.remove('../tex_papers/bloom_favar/table5.tex')
f1 = open('../tex_papers/bloom_favar/table5.tex','w')
f1.write(tmpstr)
f1.close()