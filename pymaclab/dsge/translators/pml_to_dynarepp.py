def translate(template_paramdic=None,focli=None,fpath=None):
    from pymaclab.modfiles.templates import mako_dynare_template
    
    # Render the template to be passed to dynare++
    tmp_dic = {}
    tmp_dic['focli'] = focli
    template_paramdic.update(tmp_dic)
    modstr = mako_dynare_template.render(**template_paramdic)
    
    # If a filepath has been passed then just write the Dynare++ modfile, but no more!
    if fpath != None:
        filo = open(fpath,'w')
        filo.write(modstr)
        filo.flush()
        filo.close()
        return
    else:
        return modstr