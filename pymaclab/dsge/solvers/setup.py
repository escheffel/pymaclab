# File setup.py
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration('solvers',parent_package,top_path)
    lapack = dict(get_info('lapack_opt'))
    config.add_extension('isolab', sources=['src/isolab.pyf', 
            'src/solab.f90', 'src/isolab.f90'], libraries=['lapack'],
            library_dirs=lapack['library_dirs'])
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
