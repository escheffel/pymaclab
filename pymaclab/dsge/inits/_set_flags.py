import subprocess

# Check if dynare++ is installed and in your path
dynare_flag = False
try:
    dynret = subprocess.call('dynare++',stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    dynare_flag = True
    print 'Dynare++ found on your system path, accessible from PyMacLab...'
except:
    print 'No Dynare++ installation detected on your system, unavailable...'
    dynare_flag = False