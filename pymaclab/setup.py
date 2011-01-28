#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('pymaclab', parent_package, top_path)
#    config.add_subpackage('hpfilter') # this has to go before dsge
    config.add_subpackage('dsge')
    config.add_subpackage('stats')
    config.add_extension('_hpfilter', sources=['hpfilter/hpfilt.f'])
    config.add_extension('isolab', sources=['src/isolab.pyf', 
            'src/solab.f90', 'src/isolab.f90'], libraries=['lapack'])
    return config

if __name__ == '__main__':
    print('This is the wrong setup.py file to run')
