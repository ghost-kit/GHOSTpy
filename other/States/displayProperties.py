#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
t96_dipole128vts = FindSource('t96_dipole128.vts')

# set active source
SetActiveSource(t96_dipole128vts)

# get active view
renderView2 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView2.ViewSize = [1085, 596]

# hide data in view
Hide(t96_dipole128vts, renderView2)

# set active source
SetActiveSource(t96_dipole128vts)

# show data in view
t96_dipole128vtsDisplay = Show(t96_dipole128vts, renderView2)
# trace defaults for the display properties.
t96_dipole128vtsDisplay.CubeAxesVisibility = 0
t96_dipole128vtsDisplay.Representation = 'Outline'
t96_dipole128vtsDisplay.AmbientColor = [0.16993972686350806, 0.16993972686350806, 0.16993972686350806]
t96_dipole128vtsDisplay.ColorArrayName = [None, '']
t96_dipole128vtsDisplay.DiffuseColor = [1.0, 1.0, 1.0]
t96_dipole128vtsDisplay.LookupTable = None
t96_dipole128vtsDisplay.MapScalars = 1
t96_dipole128vtsDisplay.InterpolateScalarsBeforeMapping = 1
t96_dipole128vtsDisplay.Opacity = 1.0
t96_dipole128vtsDisplay.PointSize = 2.0
t96_dipole128vtsDisplay.LineWidth = 1.0
t96_dipole128vtsDisplay.Interpolation = 'Gouraud'
t96_dipole128vtsDisplay.Specular = 0.0
t96_dipole128vtsDisplay.SpecularColor = [1.0, 1.0, 1.0]
t96_dipole128vtsDisplay.SpecularPower = 100.0
t96_dipole128vtsDisplay.Ambient = 0.0
t96_dipole128vtsDisplay.Diffuse = 1.0
t96_dipole128vtsDisplay.EdgeColor = [0.0, 0.0, 0.5]
t96_dipole128vtsDisplay.BackfaceRepresentation = 'Follow Frontface'
t96_dipole128vtsDisplay.BackfaceAmbientColor = [1.0, 1.0, 1.0]
t96_dipole128vtsDisplay.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
t96_dipole128vtsDisplay.BackfaceOpacity = 1.0
t96_dipole128vtsDisplay.Position = [0.0, 0.0, 0.0]
t96_dipole128vtsDisplay.Scale = [1.0, 1.0, 1.0]
t96_dipole128vtsDisplay.Orientation = [0.0, 0.0, 0.0]
t96_dipole128vtsDisplay.Origin = [0.0, 0.0, 0.0]
t96_dipole128vtsDisplay.Pickable = 1
t96_dipole128vtsDisplay.Texture = None
t96_dipole128vtsDisplay.Triangulate = 0
t96_dipole128vtsDisplay.NonlinearSubdivisionLevel = 1
t96_dipole128vtsDisplay.GlyphType = 'Arrow'
t96_dipole128vtsDisplay.CubeAxesColor = [0.16993972686350806, 0.16993972686350806, 0.16993972686350806]
t96_dipole128vtsDisplay.CubeAxesCornerOffset = 0.0
t96_dipole128vtsDisplay.CubeAxesFlyMode = 'Closest Triad'
t96_dipole128vtsDisplay.CubeAxesInertia = 1
t96_dipole128vtsDisplay.CubeAxesTickLocation = 'Inside'
t96_dipole128vtsDisplay.CubeAxesXAxisMinorTickVisibility = 1
t96_dipole128vtsDisplay.CubeAxesXAxisTickVisibility = 1
t96_dipole128vtsDisplay.CubeAxesXAxisVisibility = 1
t96_dipole128vtsDisplay.CubeAxesXGridLines = 0
t96_dipole128vtsDisplay.CubeAxesXTitle = 'X-Axis'
t96_dipole128vtsDisplay.CubeAxesUseDefaultXTitle = 1
t96_dipole128vtsDisplay.CubeAxesYAxisMinorTickVisibility = 1
t96_dipole128vtsDisplay.CubeAxesYAxisTickVisibility = 1
t96_dipole128vtsDisplay.CubeAxesYAxisVisibility = 1
t96_dipole128vtsDisplay.CubeAxesYGridLines = 0
t96_dipole128vtsDisplay.CubeAxesYTitle = 'Y-Axis'
t96_dipole128vtsDisplay.CubeAxesUseDefaultYTitle = 1
t96_dipole128vtsDisplay.CubeAxesZAxisMinorTickVisibility = 1
t96_dipole128vtsDisplay.CubeAxesZAxisTickVisibility = 1
t96_dipole128vtsDisplay.CubeAxesZAxisVisibility = 1
t96_dipole128vtsDisplay.CubeAxesZGridLines = 0
t96_dipole128vtsDisplay.CubeAxesZTitle = 'Z-Axis'
t96_dipole128vtsDisplay.CubeAxesUseDefaultZTitle = 1
t96_dipole128vtsDisplay.CubeAxesGridLineLocation = 'All Faces'
t96_dipole128vtsDisplay.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
t96_dipole128vtsDisplay.CustomBoundsActive = [0, 0, 0]
t96_dipole128vtsDisplay.OriginalBoundsRangeActive = [0, 0, 0]
t96_dipole128vtsDisplay.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
t96_dipole128vtsDisplay.CustomRangeActive = [0, 0, 0]
t96_dipole128vtsDisplay.UseAxesOrigin = 0
t96_dipole128vtsDisplay.AxesOrigin = [0.0, 0.0, 0.0]
t96_dipole128vtsDisplay.CubeAxesXLabelFormat = '%-#6.3g'
t96_dipole128vtsDisplay.CubeAxesYLabelFormat = '%-#6.3g'
t96_dipole128vtsDisplay.CubeAxesZLabelFormat = '%-#6.3g'
t96_dipole128vtsDisplay.StickyAxes = 0
t96_dipole128vtsDisplay.CenterStickyAxes = 0
t96_dipole128vtsDisplay.SelectionCellLabelBold = 0
t96_dipole128vtsDisplay.SelectionCellLabelColor = [0.0, 1.0, 0.0]
t96_dipole128vtsDisplay.SelectionCellLabelFontFamily = 'Arial'
t96_dipole128vtsDisplay.SelectionCellLabelFontSize = 18
t96_dipole128vtsDisplay.SelectionCellLabelItalic = 0
t96_dipole128vtsDisplay.SelectionCellLabelJustification = 'Left'
t96_dipole128vtsDisplay.SelectionCellLabelOpacity = 1.0
t96_dipole128vtsDisplay.SelectionCellLabelShadow = 0
t96_dipole128vtsDisplay.SelectionPointLabelBold = 0
t96_dipole128vtsDisplay.SelectionPointLabelColor = [1.0, 1.0, 0.0]
t96_dipole128vtsDisplay.SelectionPointLabelFontFamily = 'Arial'
t96_dipole128vtsDisplay.SelectionPointLabelFontSize = 18
t96_dipole128vtsDisplay.SelectionPointLabelItalic = 0
t96_dipole128vtsDisplay.SelectionPointLabelJustification = 'Left'
t96_dipole128vtsDisplay.SelectionPointLabelOpacity = 1.0
t96_dipole128vtsDisplay.SelectionPointLabelShadow = 0
t96_dipole128vtsDisplay.ScalarOpacityUnitDistance = 0.6681323916723285
t96_dipole128vtsDisplay.SelectMapper = 'Projected tetra'
t96_dipole128vtsDisplay.UseDataParititions = 1

# init the 'Arrow' selected for 'GlyphType'
t96_dipole128vtsDisplay.GlyphType.TipResolution = 6
t96_dipole128vtsDisplay.GlyphType.TipRadius = 0.1
t96_dipole128vtsDisplay.GlyphType.TipLength = 0.35
t96_dipole128vtsDisplay.GlyphType.ShaftResolution = 6
t96_dipole128vtsDisplay.GlyphType.ShaftRadius = 0.03
t96_dipole128vtsDisplay.GlyphType.Invert = 0

# find source
streamTracer2 = FindSource('StreamTracer2')

# hide data in view
Hide(streamTracer2, renderView2)

# set active source
SetActiveSource(streamTracer2)

# show data in view
streamTracer2Display = Show(streamTracer2, renderView2)
# trace defaults for the display properties.
streamTracer2Display.CubeAxesVisibility = 0
streamTracer2Display.Representation = 'Surface'
streamTracer2Display.AmbientColor = [0.666666666666667, 0.0817425803006027, 0.0655680170901045]
streamTracer2Display.ColorArrayName = [None, '']
streamTracer2Display.DiffuseColor = [0.0, 0.0, 1.0]
streamTracer2Display.LookupTable = None
streamTracer2Display.MapScalars = 1
streamTracer2Display.InterpolateScalarsBeforeMapping = 1
streamTracer2Display.Opacity = 1.0
streamTracer2Display.PointSize = 2.0
streamTracer2Display.LineWidth = 1.0
streamTracer2Display.Interpolation = 'Gouraud'
streamTracer2Display.Specular = 0.0
streamTracer2Display.SpecularColor = [1.0, 1.0, 1.0]
streamTracer2Display.SpecularPower = 100.0
streamTracer2Display.Ambient = 0.0
streamTracer2Display.Diffuse = 1.0
streamTracer2Display.EdgeColor = [0.0, 0.0, 0.5]
streamTracer2Display.BackfaceRepresentation = 'Follow Frontface'
streamTracer2Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
streamTracer2Display.BackfaceDiffuseColor = [1.0, 0.23099107347219, 0.267490653849088]
streamTracer2Display.BackfaceOpacity = 1.0
streamTracer2Display.Position = [0.0, 0.0, 0.0]
streamTracer2Display.Scale = [1.0, 1.0, 1.0]
streamTracer2Display.Orientation = [0.0, 0.0, 0.0]
streamTracer2Display.Origin = [0.0, 0.0, 0.0]
streamTracer2Display.Pickable = 1
streamTracer2Display.Texture = None
streamTracer2Display.Triangulate = 0
streamTracer2Display.NonlinearSubdivisionLevel = 1
streamTracer2Display.GlyphType = 'Arrow'
streamTracer2Display.CubeAxesColor = [0.666666666666667, 0.0817425803006027, 0.0655680170901045]
streamTracer2Display.CubeAxesCornerOffset = 0.0
streamTracer2Display.CubeAxesFlyMode = 'Closest Triad'
streamTracer2Display.CubeAxesInertia = 1
streamTracer2Display.CubeAxesTickLocation = 'Inside'
streamTracer2Display.CubeAxesXAxisMinorTickVisibility = 1
streamTracer2Display.CubeAxesXAxisTickVisibility = 1
streamTracer2Display.CubeAxesXAxisVisibility = 1
streamTracer2Display.CubeAxesXGridLines = 0
streamTracer2Display.CubeAxesXTitle = 'X-Axis'
streamTracer2Display.CubeAxesUseDefaultXTitle = 1
streamTracer2Display.CubeAxesYAxisMinorTickVisibility = 1
streamTracer2Display.CubeAxesYAxisTickVisibility = 1
streamTracer2Display.CubeAxesYAxisVisibility = 1
streamTracer2Display.CubeAxesYGridLines = 0
streamTracer2Display.CubeAxesYTitle = 'Y-Axis'
streamTracer2Display.CubeAxesUseDefaultYTitle = 1
streamTracer2Display.CubeAxesZAxisMinorTickVisibility = 1
streamTracer2Display.CubeAxesZAxisTickVisibility = 1
streamTracer2Display.CubeAxesZAxisVisibility = 1
streamTracer2Display.CubeAxesZGridLines = 0
streamTracer2Display.CubeAxesZTitle = 'Z-Axis'
streamTracer2Display.CubeAxesUseDefaultZTitle = 1
streamTracer2Display.CubeAxesGridLineLocation = 'All Faces'
streamTracer2Display.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
streamTracer2Display.CustomBoundsActive = [0, 0, 0]
streamTracer2Display.OriginalBoundsRangeActive = [0, 0, 0]
streamTracer2Display.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
streamTracer2Display.CustomRangeActive = [0, 0, 0]
streamTracer2Display.UseAxesOrigin = 0
streamTracer2Display.AxesOrigin = [0.0, 0.0, 0.0]
streamTracer2Display.CubeAxesXLabelFormat = '%-#6.3g'
streamTracer2Display.CubeAxesYLabelFormat = '%-#6.3g'
streamTracer2Display.CubeAxesZLabelFormat = '%-#6.3g'
streamTracer2Display.StickyAxes = 0
streamTracer2Display.CenterStickyAxes = 0
streamTracer2Display.SelectionCellLabelBold = 0
streamTracer2Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
streamTracer2Display.SelectionCellLabelFontFamily = 'Arial'
streamTracer2Display.SelectionCellLabelFontSize = 18
streamTracer2Display.SelectionCellLabelItalic = 0
streamTracer2Display.SelectionCellLabelJustification = 'Left'
streamTracer2Display.SelectionCellLabelOpacity = 1.0
streamTracer2Display.SelectionCellLabelShadow = 0
streamTracer2Display.SelectionPointLabelBold = 0
streamTracer2Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
streamTracer2Display.SelectionPointLabelFontFamily = 'Arial'
streamTracer2Display.SelectionPointLabelFontSize = 18
streamTracer2Display.SelectionPointLabelItalic = 0
streamTracer2Display.SelectionPointLabelJustification = 'Left'
streamTracer2Display.SelectionPointLabelOpacity = 1.0
streamTracer2Display.SelectionPointLabelShadow = 0
streamTracer2Display.GaussianRadius = 0.0
streamTracer2Display.ShaderPreset = 'Sphere'
streamTracer2Display.Emissive = 0
streamTracer2Display.ScaleByArray = 0
streamTracer2Display.SetScaleArray = ['POINTS', 'AngularVelocity']
streamTracer2Display.ScaleTransferFunction = 'PiecewiseFunction'
streamTracer2Display.OpacityByArray = 0
streamTracer2Display.OpacityArray = ['POINTS', 'AngularVelocity']
streamTracer2Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'Arrow' selected for 'GlyphType'
streamTracer2Display.GlyphType.TipResolution = 6
streamTracer2Display.GlyphType.TipRadius = 0.1
streamTracer2Display.GlyphType.TipLength = 0.35
streamTracer2Display.GlyphType.ShaftResolution = 6
streamTracer2Display.GlyphType.ShaftRadius = 0.03
streamTracer2Display.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
streamTracer2Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
streamTracer2Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# find source
calculator1 = FindSource('Calculator1')

# hide data in view
Hide(calculator1, renderView2)

# set active source
SetActiveSource(calculator1)

# show data in view
calculator1Display = Show(calculator1, renderView2)
# trace defaults for the display properties.
calculator1Display.CubeAxesVisibility = 0
calculator1Display.Representation = 'Outline'
calculator1Display.AmbientColor = [0.16993972686350806, 0.16993972686350806, 0.16993972686350806]
calculator1Display.ColorArrayName = ['POINTS', '']
calculator1Display.DiffuseColor = [1.0, 1.0, 1.0]
calculator1Display.LookupTable = None
calculator1Display.MapScalars = 1
calculator1Display.InterpolateScalarsBeforeMapping = 1
calculator1Display.Opacity = 1.0
calculator1Display.PointSize = 2.0
calculator1Display.LineWidth = 1.0
calculator1Display.Interpolation = 'Gouraud'
calculator1Display.Specular = 0.0
calculator1Display.SpecularColor = [1.0, 1.0, 1.0]
calculator1Display.SpecularPower = 100.0
calculator1Display.Ambient = 0.0
calculator1Display.Diffuse = 1.0
calculator1Display.EdgeColor = [0.0, 0.0, 0.5]
calculator1Display.BackfaceRepresentation = 'Follow Frontface'
calculator1Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
calculator1Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
calculator1Display.BackfaceOpacity = 1.0
calculator1Display.Position = [0.0, 0.0, 0.0]
calculator1Display.Scale = [1.0, 1.0, 1.0]
calculator1Display.Orientation = [0.0, 0.0, 0.0]
calculator1Display.Origin = [0.0, 0.0, 0.0]
calculator1Display.Pickable = 1
calculator1Display.Texture = None
calculator1Display.Triangulate = 0
calculator1Display.NonlinearSubdivisionLevel = 1
calculator1Display.GlyphType = 'Arrow'
calculator1Display.CubeAxesColor = [0.16993972686350806, 0.16993972686350806, 0.16993972686350806]
calculator1Display.CubeAxesCornerOffset = 0.0
calculator1Display.CubeAxesFlyMode = 'Closest Triad'
calculator1Display.CubeAxesInertia = 1
calculator1Display.CubeAxesTickLocation = 'Inside'
calculator1Display.CubeAxesXAxisMinorTickVisibility = 1
calculator1Display.CubeAxesXAxisTickVisibility = 1
calculator1Display.CubeAxesXAxisVisibility = 1
calculator1Display.CubeAxesXGridLines = 0
calculator1Display.CubeAxesXTitle = 'X-Axis'
calculator1Display.CubeAxesUseDefaultXTitle = 1
calculator1Display.CubeAxesYAxisMinorTickVisibility = 1
calculator1Display.CubeAxesYAxisTickVisibility = 1
calculator1Display.CubeAxesYAxisVisibility = 1
calculator1Display.CubeAxesYGridLines = 0
calculator1Display.CubeAxesYTitle = 'Y-Axis'
calculator1Display.CubeAxesUseDefaultYTitle = 1
calculator1Display.CubeAxesZAxisMinorTickVisibility = 1
calculator1Display.CubeAxesZAxisTickVisibility = 1
calculator1Display.CubeAxesZAxisVisibility = 1
calculator1Display.CubeAxesZGridLines = 0
calculator1Display.CubeAxesZTitle = 'Z-Axis'
calculator1Display.CubeAxesUseDefaultZTitle = 1
calculator1Display.CubeAxesGridLineLocation = 'All Faces'
calculator1Display.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
calculator1Display.CustomBoundsActive = [0, 0, 0]
calculator1Display.OriginalBoundsRangeActive = [0, 0, 0]
calculator1Display.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
calculator1Display.CustomRangeActive = [0, 0, 0]
calculator1Display.UseAxesOrigin = 0
calculator1Display.AxesOrigin = [0.0, 0.0, 0.0]
calculator1Display.CubeAxesXLabelFormat = '%-#6.3g'
calculator1Display.CubeAxesYLabelFormat = '%-#6.3g'
calculator1Display.CubeAxesZLabelFormat = '%-#6.3g'
calculator1Display.StickyAxes = 0
calculator1Display.CenterStickyAxes = 0
calculator1Display.SelectionCellLabelBold = 0
calculator1Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
calculator1Display.SelectionCellLabelFontFamily = 'Arial'
calculator1Display.SelectionCellLabelFontSize = 18
calculator1Display.SelectionCellLabelItalic = 0
calculator1Display.SelectionCellLabelJustification = 'Left'
calculator1Display.SelectionCellLabelOpacity = 1.0
calculator1Display.SelectionCellLabelShadow = 0
calculator1Display.SelectionPointLabelBold = 0
calculator1Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
calculator1Display.SelectionPointLabelFontFamily = 'Arial'
calculator1Display.SelectionPointLabelFontSize = 18
calculator1Display.SelectionPointLabelItalic = 0
calculator1Display.SelectionPointLabelJustification = 'Left'
calculator1Display.SelectionPointLabelOpacity = 1.0
calculator1Display.SelectionPointLabelShadow = 0
calculator1Display.ScalarOpacityUnitDistance = 0.6681323916723285
calculator1Display.SelectMapper = 'Projected tetra'
calculator1Display.UseDataParititions = 1

# init the 'Arrow' selected for 'GlyphType'
calculator1Display.GlyphType.TipResolution = 6
calculator1Display.GlyphType.TipRadius = 0.1
calculator1Display.GlyphType.TipLength = 0.35
calculator1Display.GlyphType.ShaftResolution = 6
calculator1Display.GlyphType.ShaftRadius = 0.03
calculator1Display.GlyphType.Invert = 0

# find source
contour1 = FindSource('Contour1')

# hide data in view
Hide(contour1, renderView2)

# set active source
SetActiveSource(contour1)

# show data in view
contour1Display = Show(contour1, renderView2)
# trace defaults for the display properties.
contour1Display.CubeAxesVisibility = 0
contour1Display.Representation = 'Surface'
contour1Display.AmbientColor = [0.16993972686350806, 0.16993972686350806, 0.16993972686350806]
contour1Display.ColorArrayName = [None, '']
contour1Display.DiffuseColor = [0.0, 0.0, 0.0]
contour1Display.LookupTable = None
contour1Display.MapScalars = 1
contour1Display.InterpolateScalarsBeforeMapping = 1
contour1Display.Opacity = 1.0
contour1Display.PointSize = 2.0
contour1Display.LineWidth = 1.0
contour1Display.Interpolation = 'Gouraud'
contour1Display.Specular = 0.0
contour1Display.SpecularColor = [1.0, 1.0, 1.0]
contour1Display.SpecularPower = 100.0
contour1Display.Ambient = 0.0
contour1Display.Diffuse = 1.0
contour1Display.EdgeColor = [0.0, 0.0, 0.5]
contour1Display.BackfaceRepresentation = 'Follow Frontface'
contour1Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
contour1Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
contour1Display.BackfaceOpacity = 1.0
contour1Display.Position = [0.0, 0.0, 0.0]
contour1Display.Scale = [1.0, 1.0, 1.0]
contour1Display.Orientation = [0.0, 0.0, 0.0]
contour1Display.Origin = [0.0, 0.0, 0.0]
contour1Display.Pickable = 1
contour1Display.Texture = None
contour1Display.Triangulate = 0
contour1Display.NonlinearSubdivisionLevel = 1
contour1Display.GlyphType = 'Arrow'
contour1Display.CubeAxesColor = [0.16993972686350806, 0.16993972686350806, 0.16993972686350806]
contour1Display.CubeAxesCornerOffset = 0.0
contour1Display.CubeAxesFlyMode = 'Closest Triad'
contour1Display.CubeAxesInertia = 1
contour1Display.CubeAxesTickLocation = 'Inside'
contour1Display.CubeAxesXAxisMinorTickVisibility = 1
contour1Display.CubeAxesXAxisTickVisibility = 1
contour1Display.CubeAxesXAxisVisibility = 1
contour1Display.CubeAxesXGridLines = 0
contour1Display.CubeAxesXTitle = 'X-Axis'
contour1Display.CubeAxesUseDefaultXTitle = 1
contour1Display.CubeAxesYAxisMinorTickVisibility = 1
contour1Display.CubeAxesYAxisTickVisibility = 1
contour1Display.CubeAxesYAxisVisibility = 1
contour1Display.CubeAxesYGridLines = 0
contour1Display.CubeAxesYTitle = 'Y-Axis'
contour1Display.CubeAxesUseDefaultYTitle = 1
contour1Display.CubeAxesZAxisMinorTickVisibility = 1
contour1Display.CubeAxesZAxisTickVisibility = 1
contour1Display.CubeAxesZAxisVisibility = 1
contour1Display.CubeAxesZGridLines = 0
contour1Display.CubeAxesZTitle = 'Z-Axis'
contour1Display.CubeAxesUseDefaultZTitle = 1
contour1Display.CubeAxesGridLineLocation = 'All Faces'
contour1Display.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
contour1Display.CustomBoundsActive = [0, 0, 0]
contour1Display.OriginalBoundsRangeActive = [0, 0, 0]
contour1Display.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
contour1Display.CustomRangeActive = [0, 0, 0]
contour1Display.UseAxesOrigin = 0
contour1Display.AxesOrigin = [0.0, 0.0, 0.0]
contour1Display.CubeAxesXLabelFormat = '%-#6.3g'
contour1Display.CubeAxesYLabelFormat = '%-#6.3g'
contour1Display.CubeAxesZLabelFormat = '%-#6.3g'
contour1Display.StickyAxes = 0
contour1Display.CenterStickyAxes = 0
contour1Display.SelectionCellLabelBold = 0
contour1Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
contour1Display.SelectionCellLabelFontFamily = 'Arial'
contour1Display.SelectionCellLabelFontSize = 18
contour1Display.SelectionCellLabelItalic = 0
contour1Display.SelectionCellLabelJustification = 'Left'
contour1Display.SelectionCellLabelOpacity = 1.0
contour1Display.SelectionCellLabelShadow = 0
contour1Display.SelectionPointLabelBold = 0
contour1Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
contour1Display.SelectionPointLabelFontFamily = 'Arial'
contour1Display.SelectionPointLabelFontSize = 18
contour1Display.SelectionPointLabelItalic = 0
contour1Display.SelectionPointLabelJustification = 'Left'
contour1Display.SelectionPointLabelOpacity = 1.0
contour1Display.SelectionPointLabelShadow = 0
contour1Display.GaussianRadius = 0.0
contour1Display.ShaderPreset = 'Sphere'
contour1Display.Emissive = 0
contour1Display.ScaleByArray = 0
contour1Display.SetScaleArray = [None, '']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityByArray = 0
contour1Display.OpacityArray = [None, '']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'Arrow' selected for 'GlyphType'
contour1Display.GlyphType.TipResolution = 6
contour1Display.GlyphType.TipRadius = 0.1
contour1Display.GlyphType.TipLength = 0.35
contour1Display.GlyphType.ShaftResolution = 6
contour1Display.GlyphType.ShaftRadius = 0.03
contour1Display.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(t96_dipole128vts, renderView2)

# set active source
SetActiveSource(t96_dipole128vts)

# show data in view
t96_dipole128vtsDisplay = Show(t96_dipole128vts, renderView2)

# change solid color
t96_dipole128vtsDisplay.AmbientColor = [0.6666666666666666, 0.0, 0.0]

# change representation type
t96_dipole128vtsDisplay.SetRepresentationType('Wireframe')

# change representation type
t96_dipole128vtsDisplay.SetRepresentationType('Outline')

# find view
renderView3 = FindViewOrCreate('RenderView3', viewtype='RenderView')
# uncomment following to set a specific view size
# renderView3.ViewSize = [1106, 596]

# set active view
SetActiveView(renderView3)

# find source
streamTracer3 = FindSource('StreamTracer3')

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=streamTracer3)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=streamTracer3)

# find source
plotOverLine3 = FindSource('PlotOverLine3')

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=plotOverLine3)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=plotOverLine3)

# find source
probeLocation1 = FindSource('ProbeLocation1')

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=probeLocation1)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=probeLocation1)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=streamTracer2)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=streamTracer2)

# find source
probeLocation2 = FindSource('ProbeLocation2')

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=probeLocation2)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=probeLocation2)

# find source
plotOverLine2 = FindSource('PlotOverLine2')

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=plotOverLine2)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=plotOverLine2)

# find source
slice1 = FindSource('Slice1')

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=slice1)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=slice1)

# set active source
SetActiveSource(streamTracer3)

# set active view
SetActiveView(renderView2)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=streamTracer3)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=streamTracer3)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=plotOverLine3)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=plotOverLine3)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=probeLocation1)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=probeLocation1)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=streamTracer2)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=streamTracer2)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=probeLocation2)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=probeLocation2)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=plotOverLine2)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=plotOverLine2)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=slice1)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=slice1)

# set active source
SetActiveSource(t96_dipole128vts)

# Properties modified on renderView2.AxesGrid
renderView2.AxesGrid.GridColor = [0.610910200656138, 0.035751888303959714, 0.04698252841992828]

# Properties modified on renderView2.AxesGrid
renderView2.AxesGrid.ShowGrid = 1

# Properties modified on renderView2.AxesGrid
renderView2.AxesGrid.ShowGrid = 0

# Properties modified on renderView2.AxesGrid
renderView2.AxesGrid.XLabelColor = [0.31747920958266573, 0.050141145952544444, 0.09761196307316701]

# Properties modified on renderView2.AxesGrid
renderView2.AxesGrid.XTitleColor = [0.31747920958266573, 0.06594949263752194, 0.07655451285572595]

#### saving camera placements for all active views

# current camera placement for renderView3
renderView3.InteractionMode = '2D'
renderView3.CameraPosition = [-20.0, -163.923048454133, 0.0]
renderView3.CameraFocalPoint = [-20.0, 0.0, 0.0]
renderView3.CameraViewUp = [0.0, 0.0, 1.0]
renderView3.CameraParallelScale = 38.7898992182514

# current camera placement for renderView2
renderView2.InteractionMode = '2D'
renderView2.CameraPosition = [-20.0, -163.923048454133, 0.0]
renderView2.CameraFocalPoint = [-20.0, 0.0, 0.0]
renderView2.CameraViewUp = [0.0, 0.0, 1.0]
renderView2.CameraParallelScale = 38.7898992182514

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).