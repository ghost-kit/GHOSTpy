#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active source.
dipole128vts = GetActiveSource()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1381, 992]

# get display properties
dipole128vtsDisplay = GetDisplayProperties(dipole128vts, view=renderView1)

# change representation type
dipole128vtsDisplay.SetRepresentationType('Surface')

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# set scalar coloring
ColorBy(dipole128vtsDisplay, ('POINTS', 'B'))

# rescale color and/or opacity maps used to include current data range
dipole128vtsDisplay.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
dipole128vtsDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'B'
bLUT = GetColorTransferFunction('B')
bLUT.LockDataRange = 0
bLUT.InterpretValuesAsCategories = 0
bLUT.ShowCategoricalColorsinDataRangeOnly = 0
bLUT.RescaleOnVisibilityChange = 0
bLUT.EnableOpacityMapping = 0
bLUT.RGBPoints = [2177324215807.2686, 0.231373, 0.298039, 0.752941, 2.2299980432360906e+18, 0.865003, 0.865003, 0.865003, 4.4599939091479654e+18, 0.705882, 0.0156863, 0.14902]
bLUT.UseLogScale = 0
bLUT.ColorSpace = 'Diverging'
bLUT.UseBelowRangeColor = 0
bLUT.BelowRangeColor = [0.0, 0.0, 0.0]
bLUT.UseAboveRangeColor = 0
bLUT.AboveRangeColor = [1.0, 1.0, 1.0]
bLUT.NanColor = [1.0, 1.0, 0.0]
bLUT.Discretize = 1
bLUT.NumberOfTableValues = 256
bLUT.ScalarRangeInitialized = 1.0
bLUT.HSVWrap = 0
bLUT.VectorComponent = 0
bLUT.VectorMode = 'Magnitude'
bLUT.AllowDuplicateScalars = 1
bLUT.Annotations = []
bLUT.ActiveAnnotatedValues = []
bLUT.IndexedColors = []

# get opacity transfer function/opacity map for 'B'
bPWF = GetOpacityTransferFunction('B')
bPWF.Points = [2177324215807.2686, 0.0, 0.5, 0.0, 4.4599939091479654e+18, 1.0, 0.5, 0.0]
bPWF.AllowDuplicateScalars = 1
bPWF.ScalarRangeInitialized = 1

# create a new 'Stream Tracer'
streamTracer1 = StreamTracer(Input=dipole128vts,
    SeedType='Point Source')
streamTracer1.Vectors = ['POINTS', 'B']
streamTracer1.InterpolatorType = 'Interpolator with Point Locator'
streamTracer1.SurfaceStreamlines = 0
streamTracer1.IntegrationDirection = 'BOTH'
streamTracer1.IntegratorType = 'Runge-Kutta 4-5'
streamTracer1.IntegrationStepUnit = 'Cell Length'
streamTracer1.InitialStepLength = 0.2
streamTracer1.MinimumStepLength = 0.01
streamTracer1.MaximumStepLength = 0.5
streamTracer1.MaximumSteps = 2000
streamTracer1.MaximumStreamlineLength = 20.0
streamTracer1.TerminalSpeed = 1e-12
streamTracer1.MaximumError = 1e-06
streamTracer1.ComputeVorticity = 1

# init the 'Point Source' selected for 'SeedType'
streamTracer1.SeedType.Center = [0.0, 0.0, 0.0]
streamTracer1.SeedType.NumberOfPoints = 100
streamTracer1.SeedType.Radius = 2.0

# set active source
SetActiveSource(dipole128vts)

# destroy streamTracer1
Delete(streamTracer1)
del streamTracer1

# create a new 'Cylinder'
cylinder1 = Cylinder()
cylinder1.Resolution = 6
cylinder1.Height = 1.0
cylinder1.Radius = 0.5
cylinder1.Center = [0.0, 0.0, 0.0]
cylinder1.Capping = 1

# Properties modified on cylinder1
cylinder1.Radius = 2.0

# show data in view
cylinder1Display = Show(cylinder1, renderView1)
# trace defaults for the display properties.
cylinder1Display.CubeAxesVisibility = 0
cylinder1Display.Representation = 'Surface'
cylinder1Display.AmbientColor = [1.0, 1.0, 1.0]
cylinder1Display.ColorArrayName = [None, '']
cylinder1Display.DiffuseColor = [1.0, 1.0, 1.0]
cylinder1Display.LookupTable = None
cylinder1Display.MapScalars = 1
cylinder1Display.InterpolateScalarsBeforeMapping = 1
cylinder1Display.Opacity = 1.0
cylinder1Display.PointSize = 2.0
cylinder1Display.LineWidth = 1.0
cylinder1Display.Interpolation = 'Gouraud'
cylinder1Display.Specular = 0.0
cylinder1Display.SpecularColor = [1.0, 1.0, 1.0]
cylinder1Display.SpecularPower = 100.0
cylinder1Display.Ambient = 0.0
cylinder1Display.Diffuse = 1.0
cylinder1Display.EdgeColor = [0.0, 0.0, 0.5]
cylinder1Display.BackfaceRepresentation = 'Follow Frontface'
cylinder1Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
cylinder1Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
cylinder1Display.BackfaceOpacity = 1.0
cylinder1Display.Position = [0.0, 0.0, 0.0]
cylinder1Display.Scale = [1.0, 1.0, 1.0]
cylinder1Display.Orientation = [0.0, 0.0, 0.0]
cylinder1Display.Origin = [0.0, 0.0, 0.0]
cylinder1Display.Pickable = 1
cylinder1Display.Texture = None
cylinder1Display.Triangulate = 0
cylinder1Display.NonlinearSubdivisionLevel = 1
cylinder1Display.GlyphType = 'Arrow'
cylinder1Display.CubeAxesColor = [1.0, 1.0, 1.0]
cylinder1Display.CubeAxesCornerOffset = 0.0
cylinder1Display.CubeAxesFlyMode = 'Closest Triad'
cylinder1Display.CubeAxesInertia = 1
cylinder1Display.CubeAxesTickLocation = 'Inside'
cylinder1Display.CubeAxesXAxisMinorTickVisibility = 1
cylinder1Display.CubeAxesXAxisTickVisibility = 1
cylinder1Display.CubeAxesXAxisVisibility = 1
cylinder1Display.CubeAxesXGridLines = 0
cylinder1Display.CubeAxesXTitle = 'X-Axis'
cylinder1Display.CubeAxesUseDefaultXTitle = 1
cylinder1Display.CubeAxesYAxisMinorTickVisibility = 1
cylinder1Display.CubeAxesYAxisTickVisibility = 1
cylinder1Display.CubeAxesYAxisVisibility = 1
cylinder1Display.CubeAxesYGridLines = 0
cylinder1Display.CubeAxesYTitle = 'Y-Axis'
cylinder1Display.CubeAxesUseDefaultYTitle = 1
cylinder1Display.CubeAxesZAxisMinorTickVisibility = 1
cylinder1Display.CubeAxesZAxisTickVisibility = 1
cylinder1Display.CubeAxesZAxisVisibility = 1
cylinder1Display.CubeAxesZGridLines = 0
cylinder1Display.CubeAxesZTitle = 'Z-Axis'
cylinder1Display.CubeAxesUseDefaultZTitle = 1
cylinder1Display.CubeAxesGridLineLocation = 'All Faces'
cylinder1Display.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
cylinder1Display.CustomBoundsActive = [0, 0, 0]
cylinder1Display.OriginalBoundsRangeActive = [0, 0, 0]
cylinder1Display.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
cylinder1Display.CustomRangeActive = [0, 0, 0]
cylinder1Display.UseAxesOrigin = 0
cylinder1Display.AxesOrigin = [0.0, 0.0, 0.0]
cylinder1Display.CubeAxesXLabelFormat = '%-#6.3g'
cylinder1Display.CubeAxesYLabelFormat = '%-#6.3g'
cylinder1Display.CubeAxesZLabelFormat = '%-#6.3g'
cylinder1Display.StickyAxes = 0
cylinder1Display.CenterStickyAxes = 0
cylinder1Display.SelectionCellLabelBold = 0
cylinder1Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
cylinder1Display.SelectionCellLabelFontFamily = 'Arial'
cylinder1Display.SelectionCellLabelFontSize = 18
cylinder1Display.SelectionCellLabelItalic = 0
cylinder1Display.SelectionCellLabelJustification = 'Left'
cylinder1Display.SelectionCellLabelOpacity = 1.0
cylinder1Display.SelectionCellLabelShadow = 0
cylinder1Display.SelectionPointLabelBold = 0
cylinder1Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
cylinder1Display.SelectionPointLabelFontFamily = 'Arial'
cylinder1Display.SelectionPointLabelFontSize = 18
cylinder1Display.SelectionPointLabelItalic = 0
cylinder1Display.SelectionPointLabelJustification = 'Left'
cylinder1Display.SelectionPointLabelOpacity = 1.0
cylinder1Display.SelectionPointLabelShadow = 0

# init the 'Arrow' selected for 'GlyphType'
cylinder1Display.GlyphType.TipResolution = 6
cylinder1Display.GlyphType.TipRadius = 0.1
cylinder1Display.GlyphType.TipLength = 0.35
cylinder1Display.GlyphType.ShaftResolution = 6
cylinder1Display.GlyphType.ShaftRadius = 0.03
cylinder1Display.GlyphType.Invert = 0

# Properties modified on cylinder1
cylinder1.Height = 2.0

# Properties modified on cylinder1
cylinder1.Resolution = 36

# Properties modified on renderView1
renderView1.CenterAxesVisibility = 1

# Properties modified on renderView1
renderView1.CenterAxesVisibility = 0

# set active source
SetActiveSource(dipole128vts)

# Properties modified on renderView1
renderView1.CenterAxesVisibility = 1

# Properties modified on renderView1
renderView1.CenterAxesVisibility = 0

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.Visibility = 1

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XTitle = 'X '
renderView1.AxesGrid.YTitle = 'Y '
renderView1.AxesGrid.ZTitle = 'Z '
renderView1.AxesGrid.XTitleFontFamily = 'Times'
renderView1.AxesGrid.YTitleFontFamily = 'Times'
renderView1.AxesGrid.ZTitleFontFamily = 'Times'
renderView1.AxesGrid.FacesToRender = 32
renderView1.AxesGrid.XLabelFontFamily = 'Times'
renderView1.AxesGrid.YLabelFontFamily = 'Times'
renderView1.AxesGrid.ZLabelFontFamily = 'Times'

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.FacesToRender = 4

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.FacesToRender = 6

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.AxesToLabel = 51

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.AxesToLabel = 19

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.FacesToRender = 20

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.FacesToRender = 22

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.FacesToRender = 6

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.FacesToRender = 34

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XLabelFontSize = 14
renderView1.AxesGrid.YLabelFontSize = 13

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.FacesToRender = 4

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.FacesToRender = 32

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.FacesToRender = 4

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.AxesToLabel = 16

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.AxesToLabel = 17

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.AxesToLabel = 24

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.FacesToRender = 5

# change representation type
dipole128vtsDisplay.SetRepresentationType('Wireframe')

# change representation type
dipole128vtsDisplay.SetRepresentationType('Outline')

# create a new 'Stream Tracer'
streamTracer1 = StreamTracer(Input=dipole128vts,
    SeedType='Point Source')
streamTracer1.Vectors = ['POINTS', 'B']
streamTracer1.InterpolatorType = 'Interpolator with Point Locator'
streamTracer1.SurfaceStreamlines = 0
streamTracer1.IntegrationDirection = 'BOTH'
streamTracer1.IntegratorType = 'Runge-Kutta 4-5'
streamTracer1.IntegrationStepUnit = 'Cell Length'
streamTracer1.InitialStepLength = 0.2
streamTracer1.MinimumStepLength = 0.01
streamTracer1.MaximumStepLength = 0.5
streamTracer1.MaximumSteps = 2000
streamTracer1.MaximumStreamlineLength = 20.0
streamTracer1.TerminalSpeed = 1e-12
streamTracer1.MaximumError = 1e-06
streamTracer1.ComputeVorticity = 1

# init the 'Point Source' selected for 'SeedType'
streamTracer1.SeedType.Center = [0.0, 0.0, 0.0]
streamTracer1.SeedType.NumberOfPoints = 100
streamTracer1.SeedType.Radius = 2.0

# set active source
SetActiveSource(dipole128vts)

# destroy streamTracer1
Delete(streamTracer1)
del streamTracer1

# create a new 'Stream Tracer With Custom Source'
streamTracerWithCustomSource1 = StreamTracerWithCustomSource(Input=dipole128vts,
    SeedSource=cylinder1)
streamTracerWithCustomSource1.Vectors = ['POINTS', 'B']
streamTracerWithCustomSource1.SurfaceStreamlines = 0
streamTracerWithCustomSource1.IntegrationDirection = 'BOTH'
streamTracerWithCustomSource1.IntegratorType = 'Runge-Kutta 4-5'
streamTracerWithCustomSource1.IntegrationStepUnit = 'Cell Length'
streamTracerWithCustomSource1.InitialStepLength = 0.2
streamTracerWithCustomSource1.MinimumStepLength = 0.01
streamTracerWithCustomSource1.MaximumStepLength = 0.5
streamTracerWithCustomSource1.MaximumSteps = 2000
streamTracerWithCustomSource1.MaximumStreamlineLength = 20.0
streamTracerWithCustomSource1.TerminalSpeed = 1e-12
streamTracerWithCustomSource1.MaximumError = 1e-06
streamTracerWithCustomSource1.ComputeVorticity = 1

# show data in view
streamTracerWithCustomSource1Display = Show(streamTracerWithCustomSource1, renderView1)
# trace defaults for the display properties.
streamTracerWithCustomSource1Display.CubeAxesVisibility = 0
streamTracerWithCustomSource1Display.Representation = 'Surface'
streamTracerWithCustomSource1Display.AmbientColor = [1.0, 1.0, 1.0]
streamTracerWithCustomSource1Display.ColorArrayName = ['POINTS', 'B']
streamTracerWithCustomSource1Display.DiffuseColor = [1.0, 1.0, 1.0]
streamTracerWithCustomSource1Display.LookupTable = bLUT
streamTracerWithCustomSource1Display.MapScalars = 1
streamTracerWithCustomSource1Display.InterpolateScalarsBeforeMapping = 1
streamTracerWithCustomSource1Display.Opacity = 1.0
streamTracerWithCustomSource1Display.PointSize = 2.0
streamTracerWithCustomSource1Display.LineWidth = 1.0
streamTracerWithCustomSource1Display.Interpolation = 'Gouraud'
streamTracerWithCustomSource1Display.Specular = 0.0
streamTracerWithCustomSource1Display.SpecularColor = [1.0, 1.0, 1.0]
streamTracerWithCustomSource1Display.SpecularPower = 100.0
streamTracerWithCustomSource1Display.Ambient = 0.0
streamTracerWithCustomSource1Display.Diffuse = 1.0
streamTracerWithCustomSource1Display.EdgeColor = [0.0, 0.0, 0.5]
streamTracerWithCustomSource1Display.BackfaceRepresentation = 'Follow Frontface'
streamTracerWithCustomSource1Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
streamTracerWithCustomSource1Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
streamTracerWithCustomSource1Display.BackfaceOpacity = 1.0
streamTracerWithCustomSource1Display.Position = [0.0, 0.0, 0.0]
streamTracerWithCustomSource1Display.Scale = [1.0, 1.0, 1.0]
streamTracerWithCustomSource1Display.Orientation = [0.0, 0.0, 0.0]
streamTracerWithCustomSource1Display.Origin = [0.0, 0.0, 0.0]
streamTracerWithCustomSource1Display.Pickable = 1

streamTracerWithCustomSource1Display.Triangulate = 0
streamTracerWithCustomSource1Display.NonlinearSubdivisionLevel = 1
streamTracerWithCustomSource1Display.GlyphType = 'Arrow'
streamTracerWithCustomSource1Display.CubeAxesColor = [1.0, 1.0, 1.0]
streamTracerWithCustomSource1Display.CubeAxesCornerOffset = 0.0
streamTracerWithCustomSource1Display.CubeAxesFlyMode = 'Closest Triad'
streamTracerWithCustomSource1Display.CubeAxesInertia = 1
streamTracerWithCustomSource1Display.CubeAxesTickLocation = 'Inside'
streamTracerWithCustomSource1Display.CubeAxesXAxisMinorTickVisibility = 1
streamTracerWithCustomSource1Display.CubeAxesXAxisTickVisibility = 1
streamTracerWithCustomSource1Display.CubeAxesXAxisVisibility = 1
streamTracerWithCustomSource1Display.CubeAxesXGridLines = 0
streamTracerWithCustomSource1Display.CubeAxesXTitle = 'X-Axis'
streamTracerWithCustomSource1Display.CubeAxesUseDefaultXTitle = 1
streamTracerWithCustomSource1Display.CubeAxesYAxisMinorTickVisibility = 1
streamTracerWithCustomSource1Display.CubeAxesYAxisTickVisibility = 1
streamTracerWithCustomSource1Display.CubeAxesYAxisVisibility = 1
streamTracerWithCustomSource1Display.CubeAxesYGridLines = 0
streamTracerWithCustomSource1Display.CubeAxesYTitle = 'Y-Axis'
streamTracerWithCustomSource1Display.CubeAxesUseDefaultYTitle = 1
streamTracerWithCustomSource1Display.CubeAxesZAxisMinorTickVisibility = 1
streamTracerWithCustomSource1Display.CubeAxesZAxisTickVisibility = 1
streamTracerWithCustomSource1Display.CubeAxesZAxisVisibility = 1
streamTracerWithCustomSource1Display.CubeAxesZGridLines = 0
streamTracerWithCustomSource1Display.CubeAxesZTitle = 'Z-Axis'
streamTracerWithCustomSource1Display.CubeAxesUseDefaultZTitle = 1
streamTracerWithCustomSource1Display.CubeAxesGridLineLocation = 'All Faces'
streamTracerWithCustomSource1Display.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
streamTracerWithCustomSource1Display.CustomBoundsActive = [0, 0, 0]
streamTracerWithCustomSource1Display.OriginalBoundsRangeActive = [0, 0, 0]
streamTracerWithCustomSource1Display.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
streamTracerWithCustomSource1Display.CustomRangeActive = [0, 0, 0]
streamTracerWithCustomSource1Display.UseAxesOrigin = 0
streamTracerWithCustomSource1Display.AxesOrigin = [0.0, 0.0, 0.0]
streamTracerWithCustomSource1Display.CubeAxesXLabelFormat = '%-#6.3g'
streamTracerWithCustomSource1Display.CubeAxesYLabelFormat = '%-#6.3g'
streamTracerWithCustomSource1Display.CubeAxesZLabelFormat = '%-#6.3g'
streamTracerWithCustomSource1Display.StickyAxes = 0
streamTracerWithCustomSource1Display.CenterStickyAxes = 0
streamTracerWithCustomSource1Display.SelectionCellLabelBold = 0
streamTracerWithCustomSource1Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
streamTracerWithCustomSource1Display.SelectionCellLabelFontFamily = 'Arial'
streamTracerWithCustomSource1Display.SelectionCellLabelFontSize = 18
streamTracerWithCustomSource1Display.SelectionCellLabelItalic = 0
streamTracerWithCustomSource1Display.SelectionCellLabelJustification = 'Left'
streamTracerWithCustomSource1Display.SelectionCellLabelOpacity = 1.0
streamTracerWithCustomSource1Display.SelectionCellLabelShadow = 0
streamTracerWithCustomSource1Display.SelectionPointLabelBold = 0
streamTracerWithCustomSource1Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
streamTracerWithCustomSource1Display.SelectionPointLabelFontFamily = 'Arial'
streamTracerWithCustomSource1Display.SelectionPointLabelFontSize = 18
streamTracerWithCustomSource1Display.SelectionPointLabelItalic = 0
streamTracerWithCustomSource1Display.SelectionPointLabelJustification = 'Left'
streamTracerWithCustomSource1Display.SelectionPointLabelOpacity = 1.0
streamTracerWithCustomSource1Display.SelectionPointLabelShadow = 0

# init the 'Arrow' selected for 'GlyphType'
streamTracerWithCustomSource1Display.GlyphType.TipResolution = 6
streamTracerWithCustomSource1Display.GlyphType.TipRadius = 0.1
streamTracerWithCustomSource1Display.GlyphType.TipLength = 0.35
streamTracerWithCustomSource1Display.GlyphType.ShaftResolution = 6
streamTracerWithCustomSource1Display.GlyphType.ShaftRadius = 0.03
streamTracerWithCustomSource1Display.GlyphType.Invert = 0

# hide data in view
Hide(cylinder1, renderView1)

# show color bar/color legend
streamTracerWithCustomSource1Display.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(cylinder1)

# Properties modified on cylinder1
cylinder1.Resolution = 229

# create a new 'Sphere'
sphere1 = Sphere()
sphere1.Center = [0.0, 0.0, 0.0]
sphere1.Radius = 0.5
sphere1.ThetaResolution = 8
sphere1.StartTheta = 0.0
sphere1.EndTheta = 360.0
sphere1.PhiResolution = 8
sphere1.StartPhi = 0.0
sphere1.EndPhi = 180.0

# Properties modified on sphere1
sphere1.Radius = 2.0
sphere1.ThetaResolution = 16
sphere1.PhiResolution = 4

# show data in view
sphere1Display = Show(sphere1, renderView1)
# trace defaults for the display properties.
sphere1Display.CubeAxesVisibility = 0
sphere1Display.Representation = 'Surface'
sphere1Display.AmbientColor = [1.0, 1.0, 1.0]
sphere1Display.ColorArrayName = [None, '']
sphere1Display.DiffuseColor = [1.0, 1.0, 1.0]
sphere1Display.LookupTable = None
sphere1Display.MapScalars = 1
sphere1Display.InterpolateScalarsBeforeMapping = 1
sphere1Display.Opacity = 1.0
sphere1Display.PointSize = 2.0
sphere1Display.LineWidth = 1.0
sphere1Display.Interpolation = 'Gouraud'
sphere1Display.Specular = 0.0
sphere1Display.SpecularColor = [1.0, 1.0, 1.0]
sphere1Display.SpecularPower = 100.0
sphere1Display.Ambient = 0.0
sphere1Display.Diffuse = 1.0
sphere1Display.EdgeColor = [0.0, 0.0, 0.5]
sphere1Display.BackfaceRepresentation = 'Follow Frontface'
sphere1Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
sphere1Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
sphere1Display.BackfaceOpacity = 1.0
sphere1Display.Position = [0.0, 0.0, 0.0]
sphere1Display.Scale = [1.0, 1.0, 1.0]
sphere1Display.Orientation = [0.0, 0.0, 0.0]
sphere1Display.Origin = [0.0, 0.0, 0.0]
sphere1Display.Pickable = 1
sphere1Display.Texture = None
sphere1Display.Triangulate = 0
sphere1Display.NonlinearSubdivisionLevel = 1
sphere1Display.GlyphType = 'Arrow'
sphere1Display.CubeAxesColor = [1.0, 1.0, 1.0]
sphere1Display.CubeAxesCornerOffset = 0.0
sphere1Display.CubeAxesFlyMode = 'Closest Triad'
sphere1Display.CubeAxesInertia = 1
sphere1Display.CubeAxesTickLocation = 'Inside'
sphere1Display.CubeAxesXAxisMinorTickVisibility = 1
sphere1Display.CubeAxesXAxisTickVisibility = 1
sphere1Display.CubeAxesXAxisVisibility = 1
sphere1Display.CubeAxesXGridLines = 0
sphere1Display.CubeAxesXTitle = 'X-Axis'
sphere1Display.CubeAxesUseDefaultXTitle = 1
sphere1Display.CubeAxesYAxisMinorTickVisibility = 1
sphere1Display.CubeAxesYAxisTickVisibility = 1
sphere1Display.CubeAxesYAxisVisibility = 1
sphere1Display.CubeAxesYGridLines = 0
sphere1Display.CubeAxesYTitle = 'Y-Axis'
sphere1Display.CubeAxesUseDefaultYTitle = 1
sphere1Display.CubeAxesZAxisMinorTickVisibility = 1
sphere1Display.CubeAxesZAxisTickVisibility = 1
sphere1Display.CubeAxesZAxisVisibility = 1
sphere1Display.CubeAxesZGridLines = 0
sphere1Display.CubeAxesZTitle = 'Z-Axis'
sphere1Display.CubeAxesUseDefaultZTitle = 1
sphere1Display.CubeAxesGridLineLocation = 'All Faces'
sphere1Display.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
sphere1Display.CustomBoundsActive = [0, 0, 0]
sphere1Display.OriginalBoundsRangeActive = [0, 0, 0]
sphere1Display.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
sphere1Display.CustomRangeActive = [0, 0, 0]
sphere1Display.UseAxesOrigin = 0
sphere1Display.AxesOrigin = [0.0, 0.0, 0.0]
sphere1Display.CubeAxesXLabelFormat = '%-#6.3g'
sphere1Display.CubeAxesYLabelFormat = '%-#6.3g'
sphere1Display.CubeAxesZLabelFormat = '%-#6.3g'
sphere1Display.StickyAxes = 0
sphere1Display.CenterStickyAxes = 0
sphere1Display.SelectionCellLabelBold = 0
sphere1Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
sphere1Display.SelectionCellLabelFontFamily = 'Arial'
sphere1Display.SelectionCellLabelFontSize = 18
sphere1Display.SelectionCellLabelItalic = 0
sphere1Display.SelectionCellLabelJustification = 'Left'
sphere1Display.SelectionCellLabelOpacity = 1.0
sphere1Display.SelectionCellLabelShadow = 0
sphere1Display.SelectionPointLabelBold = 0
sphere1Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
sphere1Display.SelectionPointLabelFontFamily = 'Arial'
sphere1Display.SelectionPointLabelFontSize = 18
sphere1Display.SelectionPointLabelItalic = 0
sphere1Display.SelectionPointLabelJustification = 'Left'
sphere1Display.SelectionPointLabelOpacity = 1.0
sphere1Display.SelectionPointLabelShadow = 0

# init the 'Arrow' selected for 'GlyphType'
sphere1Display.GlyphType.TipResolution = 6
sphere1Display.GlyphType.TipRadius = 0.1
sphere1Display.GlyphType.TipLength = 0.35
sphere1Display.GlyphType.ShaftResolution = 6
sphere1Display.GlyphType.ShaftRadius = 0.03
sphere1Display.GlyphType.Invert = 0

# Properties modified on sphere1
sphere1.ThetaResolution = 256

# Properties modified on sphere1
sphere1.PhiResolution = 8

# Properties modified on sphere1
sphere1.PhiResolution = 5

# set active source
SetActiveSource(streamTracerWithCustomSource1)

# Properties modified on streamTracerWithCustomSource1
streamTracerWithCustomSource1.SeedSource = sphere1

# Rescale transfer function
bLUT.RescaleTransferFunction(588.255107302, 4.45999390915e+18)

# Rescale transfer function
bPWF.RescaleTransferFunction(588.255107302, 4.45999390915e+18)

# hide data in view
Hide(sphere1, renderView1)

# set active source
SetActiveSource(sphere1)

# Properties modified on sphere1
sphere1.PhiResolution = 16

# Properties modified on sphere1
sphere1.ThetaResolution = 90

# Properties modified on sphere1
sphere1.PhiResolution = 8

# Properties modified on sphere1
sphere1.PhiResolution = 9

# Properties modified on sphere1
sphere1.PhiResolution = 25

# Properties modified on sphere1
sphere1.ThetaResolution = 64
sphere1.PhiResolution = 12

# set active source
SetActiveSource(streamTracerWithCustomSource1)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
bLUT.ApplyPreset('Black, Orange and White', True)

# convert to log space
bLUT.MapControlPointsToLogSpace()

# Properties modified on bLUT
bLUT.UseLogScale = 1

# invert the transfer function
bLUT.InvertTransferFunction()

# invert the transfer function
bLUT.InvertTransferFunction()

# invert the transfer function
bLUT.InvertTransferFunction()

# rescale color and/or opacity maps used to exactly fit the current data range
streamTracerWithCustomSource1Display.RescaleTransferFunctionToDataRange(False)

# hide color bar/color legend
streamTracerWithCustomSource1Display.SetScalarBarVisibility(renderView1, False)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
bLUT.ApplyPreset('2hot', True)

# invert the transfer function
bLUT.InvertTransferFunction()

# Properties modified on bLUT
bLUT.NumberOfTableValues = 293

# Properties modified on bLUT
bLUT.NumberOfTableValues = 1024

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [70.36400605670231, -11.558698858447174, 18.802968507360408]
renderView1.CameraFocalPoint = [-7.323403508811539, -1.1248686192780801, -1.5124514088969465]
renderView1.CameraViewUp = [-0.23488946142457565, 0.12735926758252405, 0.9636423391863373]
renderView1.CameraParallelScale = 17.320508075688775

#### uncomment the following to render all views
RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).