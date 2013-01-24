'''
.. module:: _dynarepp_parser
   :platform: Linux
   :synopsis: This is the (private) module responsible for carefully extracting information from a dynare++
   model file. This should then internally be converted to a pymaclab-style model format file.

.. moduleauthor:: Eric M. Scheffel <eric.scheffel@nottingham.edu.cn>


'''
import re
from copy import deepcopy
import numpy.matlib as mat
import numpy as np
import sys
import os
import sympycore as SP