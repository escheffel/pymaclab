# This file provides preprocess function and macros used by the
# mk_*.py scripts.
#
# Created by Pearu Peterson in Febuary 2008
#

def preprocess(source, globals_dict, tmp_cache=[1]):
    """ Preprocess source using macros from macros.py file and ones in
    globals_dict dictionary (the caller usually uses globals_dict=globals()).

    Macros start with `@` symbol and take 0 or more arguments separated
    with `;` symbol. Macros may use temporary symbols TMP, TMP0, ..., TMP5,
    that names are unique for every macro body.
    """
    g_dict = dict(globals())
    g_dict.update(globals_dict)
    result = []
    for line in source.splitlines():
        if line.lstrip().startswith('@'):
            prefix, rest = line.split('@',1)
            i = rest.index('(')
            name = rest[:i]
            tmp_cache[0] += 1
            d = {'TMP':'_tmp%s' % (tmp_cache[0])}
            for j in range(6):
                tmp_cache[0] += 1
                d['TMP'+str(j)] = '_tmp%s' % (tmp_cache[0])
            try:
                for arg in rest.strip()[i+1:-1].split(';'):
                    key, value = arg.split('=',1)
                    d[key.strip()] = value.strip()
            except Exception, msg:
                print '%s (while processing %r)' % (msg, line.lstrip())
                continue
            try:
                templ = eval(name, g_dict, {})
            except NameError, msg:
                templ = '@' + rest
                print 'NameError: %s (while processing %r)' % (msg, line.strip())
            else:
                if '@' in templ:
                    templ = preprocess(templ, globals_dict)
            result.append(prefix + '#' + rest)
            try:
                templ_d = templ % d
            except Exception, msg:
                print '%s (while processing %r)' % (msg, line.lstrip())
                #print d, `templ`
                continue
            for l in templ_d.splitlines():
                result.append(prefix + l)
        else:
            result.append(line)
    return '\n'.join(result)

NEWINSTANCE = '''\
%(OBJ)s = new(cls)
%(OBJ)s.head = %(HEAD)s
%(OBJ)s.data = %(DATA)s
'''
NEWINSTANCE = '''\
%(OBJ)s = new(cls, %(HEAD)s, %(DATA)s)
'''
NEWINSTANCE = '''\
%(OBJ)s = cls(%(HEAD)s, %(DATA)s)
'''
RETURN_NEW = '''\
@NEWINSTANCE(OBJ=%(TMP)s; HEAD=%(HEAD)s; DATA=%(DATA)s)
return %(TMP)s
'''

IF_CHECK_INT = 'if %(T)s is int or %(T)s is long:'
ELIF_CHECK_INT = 'elif %(T)s is int or %(T)s is long:'
IF_CHECK_REAL = 'if %(T)s is int or %(T)s is long or %(T)s is mpq or %(T)s is float or %(T)s is mpf:'
IF_CHECK_COMPLEX = 'if %(T)s is mpqc or %(T)s is complex or %(T)s is mpc:'

ELIF_CHECK_NUMBER = 'elif %(T)s is int or %(T)s is long or %(T)s is mpq or %(T)s is float or %(T)s is mpf or %(T)s is mpqc or %(T)s is mpc or %(T)s is complex:'

IF_CHECK_INT = 'if %(T)s in inttypes_set:'
ELIF_CHECK_INT = 'elif %(T)s in inttypes_set:'
IF_CHECK_REAL = 'if %(T)s in realtypes_set:'
IF_CHECK_COMPLEX = 'if %(T)s in complextypes_set:'
ELIF_CHECK_NUMBER = 'elif %(T)s in numbertypes_set:'

IF_CHECK_FLOAT = 'if %(T)s is float:'

ADD_TERM_VALUE_DICT='''\
%(TMP)s = %(DICT_GET)s(%(TERM)s)
if %(TMP)s is None:
    %(DICT)s[%(TERM)s] = %(USIGN)s %(VALUE)s
else:
    %(TMP)s = %(TMP)s %(SIGN)s %(VALUE)s
    if %(TMP)s:
        %(DICT)s[%(TERM)s] = %(TMP)s
    else:
        del %(DICT)s[%(TERM)s]
'''

MUL_FACTOR_VALUE_DICT='''\
%(TMP)s = %(DICT_GET)s(%(FACTOR)s)
if %(TMP)s is None:
    %(DICT)s[%(FACTOR)s] = %(USIGN)s %(VALUE)s
else:
    %(TMP)s = %(TMP)s %(SIGN)s %(VALUE)s
    if type(%(TMP)s) is cls and %(TMP)s.head is NUMBER:
        %(TMP)s = %(TMP)s.data
    if %(TMP)s:
        if %(FACTOR)s.head is NUMBER:
            del %(DICT)s[%(FACTOR)s]
            z, sym = try_power(%(FACTOR)s.data, %(TMP)s)
            if sym:
                for t1, c1 in sym:
                    @NEWINSTANCE(OBJ=tt; HEAD=NUMBER; DATA=t1)
                    @ADD_TERM_VALUE_DICT(DICT=%(DICT)s; DICT_GET=%(DICT_GET)s; TERM=tt; VALUE=c1; SIGN=+; USIGN=)
            %(NUMBER)s = %(NUMBER)s * z
        else:
            %(DICT)s[%(FACTOR)s] = %(TMP)s
    else:
        del %(DICT)s[%(FACTOR)s]
'''

