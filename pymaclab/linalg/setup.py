# File setup.py
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration('linalg',parent_package,top_path)
    config.add_data_dir('tests')
    config.add_data_dir('tests/results')
    # Let's not compile anything as this complicates things on Macintosh's
    # ordqz is still not available in Scipy for now, only qz from >=0.11
    '''
    lapack = dict(get_info('lapack_opt'))
    config.add_extension('ordqz',
                         sources = ['src/ordqz.pyf', 
                                'src/dtgsen.f'],
                         libraries=['lapack'], 
                         library_dirs=lapack['library_dirs'])
    '''
    # The code below is never needed again as qz now comes with Scipy 0.11
    '''
    config.add_extension('qz', sources = ['src/qz.pyf', 'src/dgges.f',
                        'src/zgges.f'],
                        libraries=['lapack'], 
                        library_dirs=lapack['library_dirs'])
    '''
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
