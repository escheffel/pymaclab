'''
.. module:: filters
   :platform: Linux
   :synopsis: Wrapper package module for all of the time series filters we make available in extracting the cycle and trend from
              (mostly) simulated data. But these filters are also imported into the macroeconometric class, like VAR and FAVAR to deal
              with real data.

.. moduleauthor:: Eric M. Scheffel <eric.scheffel@nottingham.edu.cn>


'''

from .hpfilter import hpfilter
from .bkfilter import bkfilter
from .cffilter import cffilter