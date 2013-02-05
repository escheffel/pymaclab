def translate(secs=None,fpath=None,dyn_vtimings={'exo':[-1,0],'endo':[-1,0],'iid':[0,1],'con':[0,1]},vardic={}):
    from pymaclab.modfiles.templates import mako_pml_template
    
    if vardic != {}:
        secs['vardic'] = vardic
    else:
        secs['vardic'] = {}
    # Render the template to be passed to dynare++
    modstr = mako_pml_template.render(**secs)