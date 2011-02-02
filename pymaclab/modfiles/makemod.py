"""
Writes a .mod file template for the user to fill in.

Notes
-----
This requires Python 2.5 or higher.
"""
import os
from string import Template

str_template = Template("""%Model Information++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$helpinfo
%Parameters+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$helpparam
%Variable Vectors+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$helpvar
%Boundary Conditions++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$helpbound
%Variable Substitution Non-Linear System++++++++++++++++++++++++++++++++++++++++
$helpsub
%Non-Linear First-Order Conditions++++++++++++++++++++++++++++++++++++++++++++++
$helpfoc
%Steady State Non-Linear System [Manual]++++++++++++++++++++++++++++++++++++++++
$helpssnl
%Steady States [Closed Form]++++++++++++++++++++++++++++++++++++++++++++++++++++
$helpsscl
%Log-Linearized Model Equations+++++++++++++++++++++++++++++++++++++++++++++++++
$helploglin
%Variance-Covariance Matrix+++++++++++++++++++++++++++++++++++++++++++++++++++++
$helpvarcov
%Minford Model Evaluation+++++++++++++++++++++++++++++++++++++++++++++++++++++++
$helpminf
%End of Model File++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
""")

nohelp=dict(zip(['helpinfo', 'helpparam', 'helpvar', 'helpbound', 
                'helpsub', 'helpfoc', 'helpssnl', 'helpsscl', 'helploglin', 
                'helpvarcov', 'helpminf'], ("",)*13))

helpinfo = """Available information: (Description, Name)
Example:
Description = This is the canonical Christiano, Eichenbaum and Evans New-Keynesian model;
Name = CEE NK Model;
"""
helpparam = """Example: 
rho = .88;
alpha_phi = .45;

Do NOT use lambda as a parameter name.
For now, make sure that you do not do integer division in the parameters. 
I.e., 1/4 will give 0 not .25.  If you want .25, do 1./4 or 1/4.
"""

def _overwrite(fname):
    overwrite = raw_input("File %s exists. Overwrite? (y/n) > " % fname) 
    overwrite = overwrite.lower()
    if overwrite != 'y' and overwrite != 'n':
        print "Please type 'y' or 'n'"
        _overwrite(fname)
    if 'y' in overwrite:
        return True
    if 'n' in overwrite:
        return False

def _get_new_name(fname):
    """
    Checks if fname exists and decide to overwrite.
    """
    overwrite = _overwrite(fname)
    if overwrite:
        return fname
    else:
        newfname = raw_input("Select a new filename (Do not respecify the path) > ")
        if not newfname.endswith('.mod'):
            newfname += '.mod'
        path = os.path.dirname(fname)
        newfname = os.path.join(path, newfname)
        if os.path.exists(newfname):
            newfname = _get_new_name(newfname)
        return newfname

def make_modfile(fname, help=False):
    """
    Creates a template .mod file for the user to edit.

    Parameters
    ----------
    fname : str
        Path and filename for new .mod file.  ".mod" extension is optional.
    help : bool, optional
        Whether or not to include some guidance within the template on how to
        edit the .mod file.
    """
    if not fname.endswith('.mod'):
        fname += ".mod"
    if os.path.exists(fname):
        fname = _get_new_name(fname)
    with open(fname, "w") as modfile:
        if not help:
            modstr = str_template.substitute(nohelp)
        else:
            raise ValueError("Help isn't written yet")
        modfile.writelines(modstr)
