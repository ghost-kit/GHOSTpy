
import paraview.simple as pv
import glob
import numpy as np

pv.LoadPlugin("/Users/jjm390/libGhostKit.dylib")

filelist=glob.glob("/Volumes/Data 1/PhD Data/MHD/LFMRCM/ToProcess/*.hdf")

filelist = sorted(filelist)

currentfile = pv.vtkLFMReader(FileNames=filelist[0])
currentfile.GridScaleFactor = 'Earth Radius: 6.5e8 cm'
currentfile.CellArrayStatus = ['Magnetic Field Vector']

calc = pv.Calculator()
calc.Input=currentfile
calc.Function = 'Magnetic Field Vector * 1e5'
calc.ResultArrayName = 'B'

DataRepresentation1 = pv.Show()
DataRepresentation1.SelectionPointFieldDataArrayName = 'Magnetic Field Vector'

segment = 16
end = None
outdir = "./outdata3/"
for file in filelist[segment:]:
    print ("Number: {0}".format(segment))
    print ("Loading File:", file)


    fileparts = str.split(file,'/')
    nameparts = str.split(fileparts[-1], '.')

    currentfile.FileNames=[file]
    # currentfile.Update()


    print ("Writing Processed File")
    # writer = pv.CreateWriter("{0}{2}_{1}.vts".format(outdir, segment, nameparts[0]), calc)
    writer = pv.CreateWriter("{0}MHD_{1}.vts".format(outdir, segment, nameparts[0]), calc)
    writer.UpdatePipeline()
    segment += 1


print ("Files Loaded")
pv.Disconnect()



# try: paraview.simple
# except: from paraview.simple import *
# paraview.simple._DisableFirstRenderCameraReset()
#
#
#
# ElkStorm_mhd_20131002T000500Z_hdf = vtkLFMReader( FileNames=['/Volumes/Data 1/PhD Data/MHD/Oct_02_2013/ElkStorm_mhd_2013-10-02T00-05-00Z.hdf'] )
#
# ElkStorm_mhd_20131002T000500Z_hdf
#
# RenderView1 = GetRenderView()
# RenderView1.CenterOfRotation = [-100495451136.0, 40808787968.0, 1104775808.0]
#
# DataRepresentation1 = Show()
# DataRepresentation1.EdgeColor = [0.0, 0.0, 0.50000762951094835]
# DataRepresentation1.SelectionPointFieldDataArrayName = 'Magnetic Field Vector'
# DataRepresentation1.ScalarOpacityUnitDistance = 3041565399.9412684
# DataRepresentation1.Representation = 'Outline'
# DataRepresentation1.ScaleFactor = 23894964428.800003
#
# RenderView1.CameraPosition = [-100495451136.0, 40808787968.0, 488924069958.68622]
# RenderView1.CameraFocalPoint = [-100495451136.0, 40808787968.0, 1104775808.0]
# RenderView1.CameraClippingRange = [480742597351.25934, 497904046961.98651]
# RenderView1.CameraParallelScale = 126256923894.66629
#
# Calculator1 = Calculator()
#
# Calculator1.Function = 'Magnetic Field Vector * 1e5'
# Calculator1.ResultArrayName = 'B'
#
# DataRepresentation2 = Show()
# DataRepresentation2.EdgeColor = [0.0, 0.0, 0.50000762951094835]
# DataRepresentation2.SelectionPointFieldDataArrayName = 'B'
# DataRepresentation2.ScalarOpacityUnitDistance = 3041565399.9412684
# DataRepresentation2.Representation = 'Outline'
# DataRepresentation2.ScaleFactor = 23894964428.800003
#
# DataRepresentation1.Visibility = 0
#
# Render()

