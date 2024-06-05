#### import the simple module from the paraview
from paraview.simple import *
import sys
import os
#include <vtkSMSettings.h>

variable2plot = 'porosity'
#variable2plot = 'composition'
#variable2plot = 'usep'
#variable2plot = 'viscosity'

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
modelpath='./'
#modelpath='../'
modelpath='/nfspm/maipe/aspect/melt_loop/outputs/'
modelname=sys.argv[1]
#modelname='rev-output-melt-k' 

print(modelpath+modelname+'/solution.pvd')
solutionpvd = PVDReader(FileName=modelpath+modelname+'/solution.pvd')
#print('whatever')
#ImportPresets(filename='/home/maipe/aspect/crust/DONE/grey-red-yellow-pink.json')
#ImportPresets(filename='/home/maipe/aspect/crust/DONE/grey-red-yellow-white.json')
ImportPresets(filename=modelpath+'/grey-rainbow.json')
ImportPresets(filename=modelpath+'/grey-red-yellow-white.json')
ImportPresets(filename=modelpath+'/grey-white-pink.json')
 
# get animation scene
animationScene1 = GetAnimationScene()
animationScene1.GoToFirst()
# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [647, 402]
renderView1.OrientationAxesVisibility = 0
renderView1.LightSwitch = 1

# show data in view
solutionpvdDisplay = Show(solutionpvd, renderView1)
# trace defaults for the display properties.
#solutionpvdDisplay.ColorArrayName = [None, '']
#solutionpvdDisplay.ScalarOpacityUnitDistance = 8662.010428739874

# reset view to fit data
renderView1.ResetCamera()

if variable2plot == 'viscosity':
  ###############################
  # PREPARE SHEAR VISCOSITY OUTPUT
  # PLOT VISCOSITY
  # set active source
  SetActiveSource(solutionpvd)
  # show data in view
  solutionpvdDisplay = Show(solutionpvd, renderView1)
  ColorBy(solutionpvdDisplay, ('POINTS', 'viscosity'))
  # rescale color and/or opacity maps used to include current data range
  solutionpvdDisplay.RescaleTransferFunctionToDataRange(False)
  # get color transfer function/color map
  viscosityLUT = GetColorTransferFunction('viscosity')
  # get opacity transfer function/opacity map
  viscosityPWF = GetOpacityTransferFunction('viscosity')
  # Rescale transfer function
  viscosityLUT.RescaleTransferFunction(1e+16, 1e+17) # THIS DOESNT WORK. WHY?
  viscosityPWF.RescaleTransferFunction(1e+16, 1e+17) # THIS DOESNT WORK. WHY?
  # convert to log space
  viscosityLUT.MapControlPointsToLogSpace()
  viscosityLUT.UseLogScale = 1
  # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
  #etaLUT.ApplyPreset('Grayscale', True)
  viscosityLUT.ApplyPreset('jet', True)
  viscosityLUT.InvertTransferFunction()
  # show/hide color bar/color legend
  solutionpvdDisplay.SetScalarBarVisibility(renderView1, False)
  # get color legend/bar for viscosityLUT in view renderView1
  viscosityLUTColorBar = GetScalarBar(viscosityLUT, renderView1)
  # Properties modified on viscosityLUTColorBar
  viscosityLUTColorBar.TitleFontSize = 10
  viscosityLUTColorBar.LabelFontSize = 10
  viscosityLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
  viscosityLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
  viscosityLUTColorBar.Position = [0.8, 0.2]

elif variable2plot == 'porosity' :
###############################
  # PREPARE POROSITY OUTPUT
  # set active source
  SetActiveSource(solutionpvd)
  # show data in view
  solutionpvdDisplay = Show(solutionpvd, renderView1)

  # set scalar coloring
  ColorBy(solutionpvdDisplay, ('POINTS', 'porosity'))
  # rescale color and/or opacity maps used to include current data range
  #solutionpvdDisplay.RescaleTransferFunctionToDataRange(True)

  # get color transfer function/color map for 'porosity'
  porosityLUT = GetColorTransferFunction('porosity')
  # get opacity transfer function/opacity map for 'porosity'
  porosityPWF = GetOpacityTransferFunction('porosity')

  # Rescale transfer function
  porosityLUT.RescaleTransferFunction(0.0, 0.5) #0.5
  # Rescale transfer function
  porosityPWF.RescaleTransferFunction(0.0, 0.5) #0.5

  # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
  porosityLUT.ApplyPreset('grey-red-yellow-white', True)
  #porosityLUT.ApplyPreset('Black-Body Radiation', True)

  # get opacity transfer function/opacity map for 'porosity'
  porosityPWF = GetOpacityTransferFunction('porosity')
  #porosityPWF.Points = [0.0, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0]
  porosityPWF.ScalarRangeInitialized = 1
  # invert the transfer function
  #porosityLUT.InvertTransferFunction()
  solutionpvdDisplay.SetScalarBarVisibility(renderView1, False)
  #get color legend/bar 
  porosityLUTColorBar = GetScalarBar(porosityLUT, renderView1)
  porosityLUTColorBar.TitleFontSize = 10
  porosityLUTColorBar.LabelFontSize = 10
  porosityLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
  porosityLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
  porosityLUTColorBar.Position = [0.8, 0.2]

elif variable2plot == 'composition':
###############################
# PREPARE COMPOSITION OUTPUT

# create a new 'Calculator'
  composition = Calculator(Input=solutionpvd)
  composition.ResultArrayName = 'composition'
  #composition.Function = '-porosity*peridotiteF+(1-porosity)*peridotite'
  composition.Function = 'porosity*peridotiteF+(1-porosity)*peridotite'

# show data in view
  compositionDisplay = Show(composition, renderView1)

# set scalar coloring
  ColorBy(compositionDisplay, ('POINTS', 'composition'))
# get color transfer function/color map
  compositionLUT = GetColorTransferFunction('composition')
  compositionPWF = GetOpacityTransferFunction('composition')
  #compositionLUT.ApplyPreset('Warm to Cool (Extended)', True)
  compositionLUT.ApplyPreset('grey-white-pink', True)
  
# rescale color and/or opacity maps used to include current data range
  #calculator1Display.RescaleTransferFunctionToDataRange(True)
  #compositionLUT.RescaleTransferFunction(-1.0, 1.0)
  compositionLUT.RescaleTransferFunction(-0.2, 2.2)
  #compositionLUT.RescaleTransferFunction(-0.2, 2.2)
  compositionLUT.InvertTransferFunction()

# show/hide color bar/color legend
  compositionDisplay.SetScalarBarVisibility(renderView1, False)

# get color legend/bar 
  cLUTColorBar = GetScalarBar(compositionLUT, renderView1)
# Properties modified 
  cLUTColorBar.TitleFontSize = 10
  cLUTColorBar.LabelFontSize = 10
  cLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
  cLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
  cLUTColorBar.Position = [0.8, 0.2]


elif variable2plot == 'usep':
###############################
# CALCULATE RELATIVE VELOCITY BETWEEN SOLID AND MELT - PREPARE OUTPUT

# create a new 'Calculator'
  usep = Calculator(Input=solutionpvd)
  usep.Function = ''
# Properties modified on calculator1
  usep.ResultArrayName = 'u_sep'
  usep.Function = '100*(u_f-velocity)'

# show data in view
  usepDisplay = Show(usep, renderView1)

# set scalar coloring
  ColorBy(usepDisplay, ('POINTS', 'u_sep'))
# get color transfer function/color map
  usepLUT = GetColorTransferFunction('u_sep')
  usepPWF = GetOpacityTransferFunction('u_sep')
  usepLUT.ApplyPreset('grey-rainbow', True)

# rescale color and/or opacity maps used to include current data range
  #usepDisplay.RescaleTransferFunctionToDataRange(True)
  usepLUT.RescaleTransferFunction(0.0, 2.0)
  #usepPWF.RescaleTransferFunction(0.0, 1.0)
  #usepLUT.InvertTransferFunction()
  usepDisplay.SetScalarBarVisibility(renderView1, False)

  # create a new 'Glyph'
  glyph1 = Glyph(Input=solutionpvd, GlyphType='2D Glyph')
  glyph1.Scalars = ['POINTS', 'T']
  #glyph1.GlyphTransform = 'Transform2'
  glyph1.Vectors = ['POINTS', 'velocity']
  # show data in view
  glyph1Display = Show(glyph1, renderView1)
  # trace defaults for the display properties.
  #glyph1Display.ColorArrayName = [None, '']
  glyph1.ScaleMode = 'vector'
  glyph1.GlyphMode = 'Every Nth Point'
  #glyph1.Stride = 400 # res5
  glyph1.Stride = 100 # res4
  glyph1.ScaleFactor = 250000.0 # 2cm/yr -> arrow 5km long
#  glyph1.ScaleFactor = 100000.0 # 5cm/yr -> arrow 5km long
  # turn off scalar coloring
  ColorBy(glyph1Display, None)



#########################################
# PREPARE TEMPERATURE OUTPUT - ISOCONTOURS
# set active source
SetActiveSource(solutionpvd)

# show data in view
solutionpvdDisplay = Show(solutionpvd, renderView1)

# reset view to fit data
#renderView1.ResetCamera()

# create a new 'Contour'
contour1 = Contour(Input=solutionpvd)
contour1.ContourBy = ['POINTS', 'T']
#contour1.Isosurfaces = [725.559814453125]
contour1.Isosurfaces = [773.0, 873.0, 973.0, 1073.0, 1173.0, 1273.0, 1373.0]
contour1.PointMergeMethod = 'Uniform Binning'
  
# Properties modified on contour1

# show data in view
contour1Display = Show(contour1, renderView1)

# Properties modified on contour1Display
contour1Display.LineWidth = 4.0

# set scalar coloring
ColorBy(contour1Display, ('POINTS', 'T'))

# get color transfer function/color map for 'T'
tLUT = GetColorTransferFunction('T')
# get opacity transfer function/opacity map for 'T'
tPWF = GetOpacityTransferFunction('T')

tLUT.ApplyPreset('Cool to Warm', True)

# hide color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1, False)

# Rescale transfer function
tLUT.RescaleTransferFunction(773.0, 1373.0)
tPWF.RescaleTransferFunction(773.0, 1373.0)


#########################################
# PREPARE P_C OUTPUT - ISOCONTOUR
# set active source
SetActiveSource(solutionpvd)

# show data in view
solutionpvdDisplay = Show(solutionpvd, renderView1)

# show p_c only where porosity>0?
# create a new 'Threshold'
thresholdPor = Threshold(Input=solutionpvd)
thresholdPor.Scalars = ['POINTS', 'porosity']
thresholdPor.ThresholdRange = [0.0, 2.0]
# create a new 'Contour'
contourPCx = Contour(Input=thresholdPor)
contourPCx.ContourBy = ['POINTS', 'p_c']
contourPCx.Isosurfaces = [-5e6]
contourPCx.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contourPC = Contour(Input=solutionpvd)
contourPC.ContourBy = ['POINTS', 'p_c']
contourPC.Isosurfaces = [-5e6]
contourPC.PointMergeMethod = 'Uniform Binning'

# show data in view
contourPCDisplay = Show(contourPC, renderView1)
contourPCDisplay.LineWidth = 4.0

# change solid color
contourPCDisplay.DiffuseColor = [0.0, 1.0, 0.0]

# hide color bar/color legend
contourPCDisplay.SetScalarBarVisibility(renderView1, False)

#############################################################
# MAKE A SNAPSHOT OR MOVIE

# set active source
if variable2plot == 'usep':
  SetActiveSource(usep)
  Hide(contour1, renderView1)
  Hide(contourPC, renderView1)
  usepDisplay.SetScalarBarVisibility(renderView1, False)
elif variable2plot == 'composition':
  SetActiveSource(composition)
  Hide(contour1, renderView1)
  Hide(contourPC, renderView1)
  #Hide(glyph1, renderView1)
  compositionDisplay.SetScalarBarVisibility(renderView1, False)
else:
  SetActiveSource(solutionpvd)
  #Hide(glyph1, renderView1)
  Hide(contour1, renderView1)
  solutionpvdDisplay.SetScalarBarVisibility(renderView1, False)

#Hide(contourPC, renderView1)


#Hide(solutionpvd, renderView1)
#Hide(calculator1, renderView1)
#Hide(composition, renderView1)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
# specifies size of plot
#renderView1.CameraPosition = [25000.0, 25000.0, 136602.5403784439]
#renderView1.CameraFocalPoint = [25000.0, 25000.0, 0.0]
#renderView1.CameraParallelScale = 26000.0 #35355.33905932738
renderView1.Background = [1, 1, 1]


# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter(Input=solutionpvd)
# show data in view
annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1)

# Properties modified on annotateTimeFilter1Display
annotateTimeFilter1Display.WindowLocation = 'AnyLocation'
annotateTimeFilter1.Scale = 1e-06
annotateTimeFilter1Display.FontSize = 35
annotateTimeFilter1.Format = '%4.1f Myr'
annotateTimeFilter1Display.Position = [0.1, 0.85]
annotateTimeFilter1Display.Position = [0.4, 0.85]
annotateTimeFilter1Display.Color = [0.0, 0.0, 0.0]

# save animation images/movie
#WriteAnimation('solution/.png', Magnification=2, FrameRate=15.0, Compression=False)
#print('here')

try: 
    os.mkdir(modelpath+'outputs/'+variable2plot+'/'+modelname+'/')
except OSError as error: 
    print(' ') 
print(modelpath+'outputs/'+variable2plot+'/'+modelname+'/'+modelname+'.png')
WriteAnimation(modelpath+'outputs/'+variable2plot+'/'+modelname+'/'+modelname+'.png', Magnification=1, FrameRate=15.0, Compression=True)

#SaveAnimation('../outputs/'+modelname+'.avi', Magnification=1, FrameRate=15.0, Compression=True) # perhaps better, yet not tested

#WriteAnimation(modelpath+'outputs/'+modelname+'.avi', Magnification=1, FrameRate=15.0, Compression=True)

renderView1.InteractionMode = '2D'
#renderView1.CameraPosition = [25000.0, 25000.0, 136602.5403784439]
#renderView1.CameraFocalPoint = [25000.0, 25000.0, 0.0]
#renderView1.CameraParallelScale = 26000.0

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
