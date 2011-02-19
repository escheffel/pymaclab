# File setup.py
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('filters',parent_package,top_path)

    config.add_extension('_hpfilter',
                         sources = ['src/hpfilter.f'])
    config.add_extension('_bkfilter',
                         sources = ['src/bkfilter.f90'])
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
