import os
import pymaclab as pm
from pymaclab import __path__ as root_dir

root_dir = root_dir[0]
fcont = os.listdir(os.path.join(root_dir,'modfiles/'))
modfiles = [x for x in fcont if x[-4:] == '.txt']
modfiles.remove('modelswloglin.txt')
modpdic = {}
for modo in modfiles:
    modpdic[modo.split('.')[0]] = os.path.join(root_dir,'modfiles',modo)
# Now expose the model files' paths to the local scope to make it accessible for users
for elem in modpdic.keys():
    exec(elem+'='+'"'+modpdic[elem]+'"')
    
# Now delete the objects we no longer want in this namespace
del os
del root_dir
del fcont
del modfiles
del modpdic
del x
del elem
del modo
del pm