#! /usr/bin/env python

import vtk, time, random, sys, os, math, signal
import mstorti.utils as utils
from numpy import *
# from scipy.io import read_array

## Read the basic configuration of the mesh. 
## Later, its node coordinates are modified
## with a vtkProgrammableFilter
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName("gascont.vtk")
reader.Update()
gridread = reader.GetOutput()

grid = vtk.vtkUnstructuredGrid()
grid.DeepCopy(gridread)
gridread = None

nnod = grid.GetNumberOfPoints()

tsleep = 0
mkvid = 0
finc = 1
Dt = 1
frame = 1
frame_end = 0
val_range = [1,-0.1]

if mkvid:
    signal.signal(signal.SIGINT,utils.my_sigint_handler)

dire = "./STEPS"
state_pattern = dire + "/gascont.state-%d.tmp.gz"
coords_pattern = dire + "/gascont-mmv.state-%d.tmp.gz"
field = 0

# Create mapper
ugmapper = vtk.vtkDataSetMapper()
ugmapper.SetInput(grid)
ugmapper.SetScalarRange(0.0,1.0)

# colorTransferFunction = vtk.vtkColorTransferFunction()
# colorTransferFunction.AddRGBPoint(0.5, 0.0, 0.0, 1.0)
# colorTransferFunction.AddRGBPoint(1.25, 0.0, 1.0, 0.0)
# colorTransferFunction.AddRGBPoint(2.0, 1.0, 0.0, 0.0)
# ugmapper.SetLookupTable(colorTransferFunction)

# Create actor
ugactor = vtk.vtkActor()
ugactor.SetMapper(ugmapper)

ugmapper2 = vtk.vtkDataSetMapper()
ugmapper2.SetInput(grid)
ugmapper2.ScalarVisibilityOff()
ugactor2 = vtk.vtkActor()
ugactor2.SetMapper(ugmapper2)
prop = ugactor2.GetProperty()
prop.SetRepresentationToWireframe()
prop.SetLineWidth(2.0)
prop.SetColor(0.0,0.0,0.0)

## ----------- CREATE THE USUAL RENDERING STUFF ---------------
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(640,480)

ren.SetBackground(.4, .4, .4)
ren.AddActor(ugactor)
ren.AddActor(ugactor2)

## ----------- WRITING TO FILE -----------------
## This is used to store the frames
## for creating a movie
w2i = vtk.vtkWindowToImageFilter()
w2i.SetInput(renWin)
w2i.Modified()
w2i.Update()

## The TIFF writer
writer = vtk.vtkTIFFWriter()
writer.SetInputConnection(w2i.GetOutputPort())
writer.SetCompressionToJPEG()
    
## ----------- MORE CAMERA SETTINGS -----------------
## Initialize camera
cam = ren.GetActiveCamera()
if None:
    cam.SetFocalPoint(0.5,0.5,1.)
    cam.SetViewUp(0.,0.,1.);
    cam.SetPosition(1.5,-0.5,3.0)
else:
    cam.SetFocalPoint(0.5,0.5,0.)
    cam.SetViewUp(0.,1.,0.);
    cam.SetPosition(0.5,0.5,2.0)

## This is added so that it gives time to set
## no border in the OpenGL window and other stuff
## like minimizing other windows. 
if mkvid:
    renWin.Render()
    raw_input("Enter something to continue: ")

## --------  THE MAIN FRAME LOOP ----------------------
# Loop while: rotating the camera and modify
# node coordinates
t = 0.0
while 1:
    t += finc*Dt
    if not mkvid and tsleep>0.0:
        time.sleep(tsleep)
    if frame % 10 == 0:
        print "frame: %d [%s]" % (frame,time.asctime())

    afile = state_pattern % frame
    dfile = utils.existz(afile)
    if dfile is None:
        print "Couldn open file %s" % afile
        sys.exit()
    ## print "processing file %s" % afile
    vals = loadtxt(dfile)
    ## vals = fromfile(dfile)
    if len(vals.shape)==2:
        nvals = vals.shape[1]
    else:
        nvals = 1
    vals = vals.reshape(nnod,nvals)

    xfile = coords_pattern % frame
    dfile = utils.existz(xfile)
    if dfile is None:
        print "Couldn open file %s" % xfile
        sys.exit()
    ## print "processing file %s" % afile
    xnod = loadtxt(dfile)
    assert len(xnod.shape)==2 and xnod.shape[1]==6

    values = grid.GetPointData().GetScalars()
    points = grid.GetPoints()
    for j in range(0,nnod):
        val = vals[j,field]
        values.SetValue(j,(val-val_range[0])/(val_range[1]-val_range[0]))
        xx = points.GetPoint(j)
        points.SetPoint(j, xnod[j,0], xnod[j,1],0*val)
    grid.Modified()
    renWin.Render()

    ## Save current frame
    if mkvid:
        assert os.path.isdir("./YUV")
        w2i.Modified()
        tiff = "./frame.tiff"
        yuv = "./YUV/frame." + str(frame) +".yuv"
        writer.SetFileName(tiff)
        writer.Write()
        os.system("convert %s %s ; gzip -f %s" % (tiff,yuv,yuv))
        os.unlink(tiff)

    if frame_end>0 and frame >= frame_end:
        break
    
    ## Update frame counter
    if mkvid:
        frame += 1
    else:
        frame += finc
