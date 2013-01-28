def translate(secs=None,fpath=None):
    from pymaclab.modfiles.templates import mako_pml_template
    
    # Render the template to be passed to dynare++
    modstr = mako_pml_template.render(**secs)