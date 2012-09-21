import os
import pymaclab as pm
from pymaclab import __path__ as root_dir
root_dir = root_dir[0]


# Now do for stable branch
def expose_stable(scopedic=None):   
	fcont = os.listdir(os.path.join(root_dir,'modfiles/models/stable'))
	modfiles = [x for x in fcont if x[-4:] == '.txt']
	modpdic = {}
	for modo in modfiles:
		modpdic[modo.split('.')[0]] = os.path.join(root_dir,'modfiles/models/stable',modo)
	# Now expose the model files' paths to the local scope to make it accessible for users
	for elem in modpdic.keys():
		exec elem+'='+'"'+modpdic[elem]+'"' in scopedic

# Now do for development branch
def expose_development(scopedic=None):
	fcont = os.listdir(os.path.join(root_dir,'modfiles/models/development'))
	modfiles = [x for x in fcont if x[-4:] == '.txt']
	modpdic = {}
	for modo in modfiles:
		modpdic[modo.split('.')[0]] = os.path.join(root_dir,'modfiles/models/development',modo)
	# Now expose the model files' paths to the local scope to make it accessible for users
	for elem in modpdic.keys():
		exec elem+'='+'"'+modpdic[elem]+'"' in scopedic

# Now do for testing branch
def expose_testing(scopedic=None):
	fcont = os.listdir(os.path.join(root_dir,'modfiles/models/testing'))
	modfiles = [x for x in fcont if x[-4:] == '.txt']
	modpdic = {}
	for modo in modfiles:
		modpdic[modo.split('.')[0]] = os.path.join(root_dir,'modfiles/models/testing',modo)
	# Now expose the model files' paths to the local scope to make it accessible for users
	for elem in modpdic.keys():
		exec elem+'='+'"'+modpdic[elem]+'"' in scopedic


# Now delete the objects we no longer want in this namespace
def delete_all():
	if 'os' in dir(): del os
	if 'root_dir' in dir(): del root_dir
	if 'fcont' in dir(): del fcont
	if 'modfiles' in dir(): del modfiles
	if 'modpdic' in dir(): del modpdic
	if 'x' in dir(): del x
	if 'elem' in dir(): del elem
	if 'modo' in dir(): del modo
	if 'pm' in dir(): del pm
