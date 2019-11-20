#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
t96_128_dst100vts = FindSource('t96_128_dst100.vts')

# create a new 'Calculator'
calculator1 = Calculator(Input=t96_128_dst100vts)
calculator1.Function = ''

# Properties modified on calculator1
calculator1.ResultArrayName = 'Bmag'
calculator1.Function = 'mag(B)'

# get active view
renderView2 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView2.ViewSize = [1085, 596]

# show data in view
calculator1Display = Show(calculator1, renderView2)
# trace defaults for the display properties.
calculator1Display.Representation = 'Outline'
calculator1Display.ColorArrayName = ['POINTS', '']
calculator1Display.GlyphType = 'Arrow'
calculator1Display.ScalarOpacityUnitDistance = 0.6681323916723285

# hide data in view
Hide(t96_128_dst100vts, renderView2)

# find source
probeLocation2 = FindSource('ProbeLocation2')

# find source
probeLocation1 = FindSource('ProbeLocation1')

# find source
plotOverLine2 = FindSource('PlotOverLine2')

# find source
plotOverLine3 = FindSource('PlotOverLine3')

# find source
streamTracer3 = FindSource('StreamTracer3')

# find source
streamTracer2 = FindSource('StreamTracer2')

# find source
t96_dipole128vts = FindSource('t96_dipole128.vts')

# hide data in view
Hide(calculator1, renderView2)

# set active source
SetActiveSource(calculator1)

# show data in view
calculator1Display = Show(calculator1, renderView2)

# hide data in view
Hide(calculator1, renderView2)

# show data in view
calculator1Display = Show(calculator1, renderView2)

# create a new 'Slice'
slice1 = Slice(Input=calculator1)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [-20.0, 0.0, 0.0]

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 1.0, 0.0]

# get color transfer function/color map for 'Bmag'
bmagLUT = GetColorTransferFunction('Bmag')

# show data in view
slice1Display = Show(slice1, renderView2)
# trace defaults for the display properties.
slice1Display.ColorArrayName = ['POINTS', 'Bmag']
slice1Display.LookupTable = bmagLUT
slice1Display.GlyphType = 'Arrow'
slice1Display.SetScaleArray = ['POINTS', 'Bmag']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['POINTS', 'Bmag']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView2, True)

# get opacity transfer function/opacity map for 'Bmag'
bmagPWF = GetOpacityTransferFunction('Bmag')

# hide data in view
Hide(slice1, renderView2)

# create a new 'Contour'
contour1 = Contour(Input=slice1)
contour1.ContourBy = ['POINTS', 'Bmag']
contour1.Isosurfaces = [3243986.8102642912]
contour1.PointMergeMethod = 'Uniform Binning'

# Properties modified on contour1
contour1.Isosurfaces = [300.0, 200.0, 100.0]

# show data in view
contour1Display = Show(contour1, renderView2)
# trace defaults for the display properties.
contour1Display.ColorArrayName = [None, '']
contour1Display.GlyphType = 'Arrow'
contour1Display.SetScaleArray = [None, '']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = [None, '']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(slice1, renderView2)

# change solid color
contour1Display.DiffuseColor = [0.0, 0.0, 0.0]

#### saving camera placements for all active views

# current camera placement for renderView2
renderView2.InteractionMode = '2D'
renderView2.CameraPosition = [-20.0, -163.923048454133, 0.0]
renderView2.CameraFocalPoint = [-20.0, 0.0, 0.0]
renderView2.CameraViewUp = [0.0, 0.0, 1.0]
renderView2.CameraParallelScale = 38.7898992182514

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).