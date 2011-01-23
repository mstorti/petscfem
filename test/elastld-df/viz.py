#! /usr/bin/env python

import vtk, time, random, sys, os, math, signal
import mstorti.utils as utils
from numpy import *
# from scipy.io import read_array

if None:
    v = vtk.vtkGraphicsFactory()
    v.SetUseMesaClasses(1)
    v.SetOffScreenOnlyMode(1)
    del v
    w = vtk.vtkImagingFactory()
    w.SetUseMesaClasses(1)
    del w

## Read the basic configuration of the mesh. 
## Later, its node coordinates are modified
## with a vtkProgrammableFilter
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName("beam.vtk")
reader.Update()
gridread = reader.GetOutput()

grid = vtk.vtkUnstructuredGrid()
grid.DeepCopy(gridread)
points0 = gridread.GetPoints()
points = grid.GetPoints()

nnod = grid.GetNumberOfPoints()

# global frame_forced, frame
utils.frame_forced = None

ndim = 3
tsleep = 0
mkvid = 1
finc = 1
state_inc = -2
Dt = 1
frame_start = 0
frame_end = 1600
valrange = [0.6,0.0]

if mkvid:
    signal.signal(signal.SIGINT,utils.my_sigint_handler)

dire = "./STEPSc"
state_pattern = dire + "/elastld.state-%d.tmp.gz"
field = 0

ugmapper3 = vtk.vtkDataSetMapper()
ugmapper3.SetInput(grid)
ugmapper3.ScalarVisibilityOff()
ugactor3 = vtk.vtkActor()
ugactor3.SetMapper(ugmapper3)
prop = ugactor3.GetProperty()
# prop.SetOpacity(0.8)

ugmapper4 = vtk.vtkDataSetMapper()
ugmapper4.SetInput(grid)
ugmapper4.ScalarVisibilityOff()
ugactor4 = vtk.vtkActor()
ugactor4.SetMapper(ugmapper4)
prop = ugactor4.GetProperty()
prop.SetRepresentationToWireframe()
prop.SetLineWidth(1.6)
prop.SetColor(0.0,0.0,0.0)

## ----------- CREATE THE USUAL RENDERING STUFF ---------------
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(640,480)
# renWin.OffScreenRenderingOn();
# renWin.SetOffScreenRendering(1);

ren.SetBackground(.4, .4, .4)
ren.AddActor(ugactor3)
ren.AddActor(ugactor4)

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
if ndim==3:
    to = array([0.,3.,0.])
    cam.SetFocalPoint(to)
    cam.SetViewUp(0.,1.,0.);
else:
    to = array([0.,0.,0.])
    cam.SetFocalPoint(to)
    cam.SetViewUp(0.,1.,0.);
    cam.SetPosition(0.,0.,2.3)

def draw_axes(xc,L,D):
    cyl = vtk.vtkCylinderSource()
    cyl.SetHeight(L)
    #cyl.SetCenter(xc)
    cyl.SetRadius(D/2.0)
    cyl.SetResolution(30)

    cylMapper = vtk.vtkPolyDataMapper()
    cylMapper.SetInputConnection(cyl.GetOutputPort() )

    ## Create the cylinder actor
    cylY = vtk.vtkActor()
    cylY.SetMapper(cylMapper)
    cylY.SetPosition(xc)

    cylX = vtk.vtkActor()
    cylX.SetMapper(cylMapper)
    cylX.RotateZ(90.0)
    cylX.SetPosition(xc)

    cylZ = vtk.vtkActor()
    cylZ.SetMapper(cylMapper)
    cylZ.RotateX(90.0)
    cylZ.SetPosition(xc)

    cylX.GetProperty().SetColor(1.0,0.0,0.0)
    cylY.GetProperty().SetColor(0.0,1.0,0.0)
    cylZ.GetProperty().SetColor(0.0,0.0,1.0)
    return [cylX,cylY,cylZ]

## Draw axes
axes = draw_axes([0.,0.,0.],30,0.1)
ren.AddActor(axes[0])
ren.AddActor(axes[1])
ren.AddActor(axes[2])

wcam = -0.005
## Initial position angle for the cam
phi0 = -75.0*pi/180.0
## The speed with which the radius varies
wRcam = 1.33*wcam
Rcammax = 10
Rcammin = 0.7*Rcammax
## Height of camera
Hcam_fac = 0.5
## height of point where the camera points to
Hto = 3.0

## This is added so that it gives time to set
## no border in the OpenGL window and other stuff
## like minimizing other windows. 
if mkvid:
    renWin.Render()
    renWin.BordersOff()
#    raw_input("Enter something to continue: ")

def load_state(afile):
    dfile = utils.existz(afile)
    global frame
    if dfile is None:
        print "Couldn open file %s" % afile
        if mkvid:
            sys.exit()
        else:
            frame = frame_start
            return None
    else:
        vals = loadtxt(dfile)
        return vals

## --------  THE MAIN FRAME LOOP ----------------------
# Loop while: rotating the camera and modify
# node coordinates
frame = frame_start
t = 0.0
vals_0 = None
vals_1 = None
while 1:
    t += finc*Dt
    if utils.frame_forced:
        frame = utils.frame_forced
        utils.frame_forced = None
    if not mkvid and tsleep>0.0:
        time.sleep(tsleep)
    if frame % 10 == 0:
        print "frame: %d [%s]" % (frame,time.asctime())

    if ndim==3:
        phi = phi0+ wcam*frame
        Rcamav =(Rcammin+Rcammax)/2.0
        DRcam = (Rcammin-Rcammax)/2.0
        Rcam = Rcamav + DRcam*sin(wRcam*frame)
        Hcam = Hcam_fac*Rcam
        view_dir = array([Rcam*cos(phi),Hcam,Rcam*sin(phi)])
        pfrom = to + view_dir
        cam.SetPosition(pfrom)

    renWin.Render()
    ren.ResetCameraClippingRange()

    if state_inc>0:
        valsnew = load_state(state_pattern % (frame*state_inc))
        if valsnew != None:
            vals = valsnew
        else:
            frame = frame_start
    else:
        if vals_0 == None:
            vals_1 = loadtxt(state_pattern % 0)
            vals_0 = vals_1
            next_state = 1
        if frame % (-state_inc)==0:
            vals_0 = vals_1
            vals_1 = load_state(state_pattern % next_state)
            next_state += 1
        alpha = float(frame % (-state_inc))/(-state_inc)
        vals = (1-alpha)*vals_0 + alpha*vals_1
        
    ## vals = fromfile(dfile)
    if len(vals.shape)==2:
        nvals = vals.shape[1]
    else:
        nvals = 1
        vals = vals.reshape(nnod,nvals)

    values = grid.GetPointData().GetScalars()
    for j in range(0,nnod):
        xx = points0.GetPoint(j)
        points.SetPoint(j, xx[0]+vals[j,0], xx[1]+vals[j,1], xx[2]+vals[j,2])
    grid.Modified()
    # renWin.Render()

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

    if frame_end>=0 and frame >= frame_end:
        break
    
    ## Update frame counter
    if mkvid:
        frame += 1
    else:
        frame += finc
