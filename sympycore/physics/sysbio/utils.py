import types
import sympycore

def objsize(obj):
    """ Recursively compute memory size of the object in bytes.

    Returns
    -------
    size : int
      Number of bytes that the object consumes in memory.
    """
    import numpy
    if obj is None:
        return 0 # because None is singleton
    if isinstance (obj, numpy.ndarray):
        return 100+obj.nbytes
    if isinstance (obj, (int, float, bool)): 
        return 24
    if isinstance (obj, (long,complex)): 
        return 32
    if isinstance (obj, tuple): 
        sz = 56
        sz += len (obj)*8
        sz += sum (map (objsize, obj))
        return sz
    if isinstance (obj, list): 
        sz = 136 if obj else 72
        sz += sum (map (objsize, obj))
        return sz
    if isinstance (obj, str):
        sz = 40
        sz += (len(obj) // 8)*8 if obj else 0
        return sz
    if isinstance (obj, dict):
        sz = 280
        for k,v in obj.iteritems ():
            sz += objsize(k) + objsize(v)
        return sz
    if isinstance (obj, set):
        sz = 232
        sz += sum (map (objsize, obj))
        return sz
    if isinstance(obj, types.InstanceType):
        sz = 72
        for k,v in obj.__dict__.iteritems ():
            sz += objsize(k) + objsize(v)
        return sz
    if isinstance(obj, object):
        sz = 64
        for k,v in obj.__dict__.iteritems ():
            sz += objsize(k) + objsize(v)
        return sz
    if obj is None:
        return 16
    raise NotImplementedError ('objsize for %s' % (type (obj)))

def obj2num(s, abs_tol=1e-16):
    if isinstance(s, str):
        f = eval(s)
    else:
        f = s
    #return float(f) # will induce numerical errors and incorrect rank for GJE algorithm.
    i = int (f)
    if i==f:
        return i
    return sympycore.f2q(f, abs_tol)

def time2str(s, last_unit='s'):
    """ Return human readable time string from seconds.

    Examples
    --------
    >>> from iocbio.utils import time_to_str
    >>> print time_to_str(123000000)
    3Y10M24d10h40m
    >>> print time_to_str(1230000)
    14d5h40m
    >>> print time_to_str(1230)
    20m30.0s
    >>> print time_to_str(0.123)
    123ms
    >>> print time_to_str(0.000123)
    123us
    >>> print time_to_str(0.000000123)
    123ns

    """
    seconds_in_year = 31556925.9747 # a standard SI year
    orig_s = s
    years = int(s / (seconds_in_year))
    r = []
    if years:
        r.append ('%sY' % (years))
        s -= years * (seconds_in_year)
    if last_unit=='Y': s = 0
    months = int(s / (seconds_in_year/12.0))
    if months:
        r.append ('%sM' % (months))
        s -= months * (seconds_in_year/12.0)
    if last_unit=='M': s = 0
    days = int(s / (60*60*24))
    if days:
        r.append ('%sd' % (days))
        s -= days * 60*60*24
    if last_unit=='d': s = 0
    hours = int(s / (60*60))
    if hours:
        r.append ('%sh' % (hours))
        s -= hours * 60*60
    if last_unit=='h': s = 0
    minutes = int(s / 60)
    if minutes:
        r.append ('%sm' % (minutes))
        s -= minutes * 60
    if last_unit=='m': s = 0
    seconds = int(s)
    if seconds:
        r.append ('%ss' % (seconds))
        s -= seconds
    if last_unit=='s': s = 0
    mseconds = int(s*1000)
    if mseconds:    
        r.append ('%sms' % (mseconds))
        s -= mseconds / 1000    
    if last_unit=='ms': s = 0
    useconds = int(s*1000000)
    if useconds:
        r.append ('%sus' % (useconds))
        s -= useconds / 1000000
    if last_unit=='us': s = 0
    nseconds = int(s*1000000000)
    if nseconds:
        r.append ('%sns' % (nseconds))
        s -= nseconds / 1000000000
    if not r:
        return '0'
    return ''.join(r)
