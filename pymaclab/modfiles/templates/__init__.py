try:
    import wheezy
    from wheezy_dynare_template import wheezy_dynare_template
    from wheezy_template import wheezy_template
except:
    pass

try:
    import jinja2
    from jinja2_template import jinja2_template
except:
    pass

try:
    import cheetah
    from cheetah_template import cheetah_template
except:
    pass
  
try:
    import mako
    from mako_dynare_template import mako_dynare_template
except:
    pass

try:
    import mako
    from mako_pml_template import mako_pml_template
except:
    pass