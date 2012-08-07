#!/usr/bin/env python

"""
Generates test files for testing of import into Maya.

USAGE:  ./MakeTestFiles.py
      
$Id: MakeTestFiles.py,v 1.2 2008/04/24 22:17:26 rsz Exp $
"""

import vtk
import sys

#######################################
# For visual inspection.

ren = vtk.vtkRenderer()
win = vtk.vtkRenderWindow()
win.AddRenderer(ren)
win.SetSize(300, 300)
win.SetPosition(500,500)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(win)
style = vtk.vtkInteractorStyleTrackballCamera()
iren.SetInteractorStyle(style)
ren.SetBackground(.4, .5, .6)

# Axes indicator.
ax = vtk.vtkAxesActor()
ax.SetShaftTypeToCylinder()
ax.SetXAxisLabelText("X")
ax.SetYAxisLabelText("Y")
ax.SetZAxisLabelText("Z")
ax.SetTotalLength(1,1,1)

ren.AddActor(ax)

#######################################
# Structured Grid -> Maya Paricles
pts = vtk.vtkPoints()
n = 10
c=0
for i in xrange(n):
    for j in xrange(n):
        for k in xrange(n):
            pts.InsertPoint(c, i,j,k)
            c += 1

scalars = vtk.vtkFloatArray()
scalars.SetNumberOfTuples(n**3)
for i in xrange(int(n**3/3.)): scalars.SetTuple1(i,.1)
for i in xrange(int(n**3/3.), int(2*(n**3)/3.)): scalars.SetTuple1(i,.5)
for i in xrange(int(2*(n**3)/3.), int(n**3)): scalars.SetTuple1(i,.9)

# Wrap as a structured grid even though the data is not.
sg = vtk.vtkStructuredGrid()
sg.GetPointData().SetScalars(scalars)
sg.SetPoints(pts)
sg.SetDimensions(n,n,n)

w = vtk.vtkDataSetWriter()
w.SetFileTypeToASCII()
w.SetFileName('tmp-sgrid.vtk')
w.SetInput(sg)
w.Write()
w.Update()

geom = vtk.vtkStructuredGridGeometryFilter()
geom.SetInput(sg)
geom.Update()

m = vtk.vtkPolyDataMapper()
m.SetInput(geom.GetOutput())

a = vtk.vtkActor()
a.SetMapper(m)
a.GetProperty().SetRepresentationToPoints()
a.GetProperty().SetPointSize(15)

import time

ren.AddActor(a)

iren.Initialize()
win.Render()
time.sleep(2)

#######################################
# ImageData -> Fluid
img = vtk.vtkImageData()
img.GetPointData().SetScalars(scalars)
img.SetSpacing(.5, .75, 1)
img.SetDimensions(n,n,n)

w.SetFileName('tmp-img.vtk')
w.SetInput(img)
w.Write()
w.Update()

geom = vtk.vtkImageDataGeometryFilter()
geom.SetInput(img)

m.SetInput(geom.GetOutput())
win.Render()
time.sleep(2)

#######################################
# PolyData -> Maya Mesh
from math import sin, sqrt

pts = vtk.vtkPoints()
scalars = vtk.vtkFloatArray()

x = 0. 
xmax = 10.
eps = 1e-6
delta = .5
while x < xmax:
    y = 0.
    while y < xmax:
        z = sqrt( (x+eps)**2 + (y+eps)**2 )
        z = 5*sin(z)/z
        pts.InsertNextPoint(x,y,z)
        scalars.InsertNextTuple1(z)
        y += delta
    x += delta

sg = vtk.vtkStructuredGrid()
sg.SetPoints(pts)
sg.GetPointData().SetScalars(scalars)
sg.SetDimensions( int(xmax/delta),
                  int(xmax/delta),
                  1)

scalars.ComputeRange(0)
lut = vtk.vtkLookupTable()
# Insert some transparency for testing.
lut.SetAlphaRange(.25, 1)
lut.ForceBuild()
lut.SetRange(scalars.GetRange())                      
scalars.SetLookupTable(lut)

geom = vtk.vtkStructuredGridGeometryFilter()
geom.SetInput(sg)
geom.Update()
    
poly = geom.GetOutput()

w.SetFileName('tmp-polyspt.vtk')
w.SetInput(poly)
w.Write()
w.Update()

mapper = vtk.vtkPolyDataMapper()
mapper.SetInput(poly)
mapper.SetLookupTable(lut)
mapper.UseLookupTableScalarRangeOn()

ren.RemoveActor(a)

a = vtk.vtkActor()
a.SetMapper(mapper)
ren.AddActor(a)
win.Render()
time.sleep(2)

#######################################
# PolyData Lines -> Maya Linear (degree 1) Lines

quadric = vtk.vtkQuadric()
quadric.SetCoefficients(.5, 1, .2, 0, .1, 0, 0, .2, 0, 0)

sample = vtk.vtkSampleFunction()
sample.SetSampleDimensions(30, 30, 30)
sample.SetImplicitFunction(quadric)
sample.ComputeNormalsOff()
sample.Update()

extract = vtk.vtkExtractVOI()
extract.SetInputConnection(sample.GetOutputPort())
extract.SetVOI(0, 29, 0, 29, 15, 15)
extract.SetSampleRate(1, 2, 3)

contours = vtk.vtkContourFilter()
contours.SetInputConnection(extract.GetOutputPort())
contours.GenerateValues(4, 0.0, 1.2)
contours.Update()

poly = contours.GetOutput()
lut.SetRange(0, 1.2)
poly.GetPointData().GetScalars().SetLookupTable(lut)

# Reduce number of necessary lines by finding connected lines.
strip = vtk.vtkStripper()
strip.SetInput(poly)
strip.Update()

poly = strip.GetOutput()
w.SetFileName('tmp-polyslines.vtk')
w.SetInput(poly)
w.Write()
w.Update()

contMapper = vtk.vtkPolyDataMapper()
contMapper.SetInputConnection(contours.GetOutputPort())
contMapper.UseLookupTableScalarRangeOn()
contMapper.SetLookupTable(lut)

ren.RemoveActor(a)
a = vtk.vtkActor()
a.SetMapper(contMapper)
a.GetProperty().SetLineWidth(3)

ren.AddActor(a)
win.Render()
time.sleep(2)

###############################################
# UV Textured PolyData -> Maya Mesh
ren.RemoveActor(a)

pts = vtk.vtkPoints()
uv = vtk.vtkFloatArray()
uv.SetNumberOfComponents(2)

for i in xrange(10):
    for j in xrange(10):
        pts.InsertNextPoint(i/9., j/9., 0)
        d,m = divmod(j, 5)
        u = m/(5.-1)
        d,m = divmod(i, 5)
        v = m/(5.-1)
        uv.InsertNextTuple2(u, v)

# Wrap as a structured grid.
sg = vtk.vtkStructuredGrid()
sg.GetPointData().SetTCoords(uv)
sg.SetPoints(pts)
sg.SetDimensions(10,10,1)

geom = vtk.vtkStructuredGridGeometryFilter()
geom.SetInput(sg)
geom.Update()

w = vtk.vtkDataSetWriter()
w.SetFileTypeToASCII()
w.SetFileName('tmp-polytextured.vtk')
w.SetInput(geom.GetOutput())
w.Write()

m = vtk.vtkPolyDataMapper()
m.SetInput(geom.GetOutput())

a = vtk.vtkActor()
a.SetMapper(m)

ren.AddActor(a)
win.Render()

###############################################
# PolyData with normals -> Maya Mesh
def GetNormalLines(vtx, normals, scale):
    """
    Returns polydata with lines to visualize normals.
    This can be used for either vertex or face normals;
    just specify the normals start points with "vtx".
    """
    linePoints = vtk.vtkPoints()
    aLine = vtk.vtkLine()
    aLineGrid = vtk.vtkUnstructuredGrid()
    aLineGrid.Allocate(1, 1)
    aLineGrid.SetPoints(linePoints)

    for i in xrange(vtx.GetNumberOfTuples()):
        xyz0 = vtx.GetTuple3(i)
        nxyz = normals.GetTuple3(i)
        linePoints.InsertNextPoint(xyz0[0], xyz0[1], xyz0[2])
        linePoints.InsertNextPoint(xyz0[0]+nxyz[0]*scale,
                                   xyz0[1]+nxyz[1]*scale,
                                   xyz0[2]+nxyz[2]*scale)
        aLine.GetPointIds().SetId(0, 2*i)
        aLine.GetPointIds().SetId(1, 2*i+1)
        aLineGrid.InsertNextCell(aLine.GetCellType(),
                                 aLine.GetPointIds())
    
    return aLineGrid

ren.RemoveActor(a)

sphere = vtk.vtkSphereSource()
sphere.Update()

normals = vtk.vtkPolyDataNormals()
normals.SetInput(sphere.GetOutput())
normals.Update()

n = normals.GetOutput().GetPointData().GetNormals()
for i in xrange(n.GetNumberOfTuples()):
    xyz = n.GetTuple3(i)
    n.SetTuple3(i, xyz[0], xyz[1], 0.)

w = vtk.vtkDataSetWriter()
w.SetFileTypeToASCII()
w.SetFileName('tmp-polynormals.vtk')
w.SetInput(normals.GetOutput())
w.Write()

m = vtk.vtkPolyDataMapper()
m.SetInput(normals.GetOutput())

a = vtk.vtkActor()
a.SetMapper(m)

aLineGrid = GetNormalLines(normals.GetOutput().GetPoints().GetData(), n, .25)

aLineMapper = vtk.vtkDataSetMapper()
aLineMapper.SetInput(aLineGrid)

aLineActor = vtk.vtkActor()
aLineActor.SetMapper(aLineMapper)
aLineActor.GetProperty().SetDiffuseColor(1, 1, 0)

ren.AddActor(a)
ren.AddActor(aLineActor)
win.Render()

###############################################
iren.Start()
