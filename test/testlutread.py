import vtk

"""
This shows that the lookup table export for scalars in a polydata
is flawed.  The scalar range doesn't agree with the LUT before file
export.
"""

r = vtk.vtkDataSetReader()
r.SetFileName('tmp-polyspt.vtk')
r.Update()

print r.GetOutput()
