"""
20080525 rsz Created and tested example.

$ Id: $
"""
import maya

def importvtk(*args, **kargs):
    """
    Only one non-keyword argument is optional, and it's the VTK Python
    object variable name.  All other input arguments must be keyword arguments.
    See the importvtk help.

    # Example
    import vtk
    c = vtk.vtkConeSource()
    c.Update()
    p = c.GetOutput()
    from importvtk import importvtk
    importvtk('p')
    importvtk(r=p.__this__)
    importvtk(reference=p.__this__)
    """
    if not maya.cmds.pluginInfo('importvtk', query=True, loaded=True):
        print 'ERROR: importvtk plug-in not loaded.  Did you compile "importvtk.mll" and place it in your MAYA_PLUG_IN_PATH?  A pure Python implementation of importvtk isn\'t available.'
        return
    
    kargstr = ''
    if len(kargs):
        for k in kargs:
            kstr = repr(kargs[k])
            if (kstr[0] == '\'') and (kstr[-1] == '\''):
                # MEL doesn't like single quoted strings, so make double.
                kstr = '"' + kstr[1:-1] + '"'
                
            kargstr += ' -%s %s ' % (k, kstr)
    
    if len(args) > 0:        
        maya.mel.eval('string $r = `python "%s.__this__"`; importvtk -r $r %s;' % (args[0], kargstr))
    else:
        maya.mel.eval('importvtk %s;' % (kargstr))
