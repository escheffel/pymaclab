def translate(template_paramdic=None,fpath=None):
    import pymaclab.modfiles.templates.wheezy_template as template
    
    # Render the template to be passed to wheezy pymaclab template
    modstr = template.render(template_paramdic)
    
    # If a filepath has been passed then just write the pymaclab modfile, but no more!
    if fpath != None:
        filo = open(fpath,'w')
        filo.write(modstr)
        filo.flush()
        filo.close()
        return
    else:
        return modstr