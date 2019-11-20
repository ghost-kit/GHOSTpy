import vtk

file = vtk.vtkXMLStructuredGridReader()
file.SetFileName("out/fields/dipole_dp0_3_small_grid.vts")
file.Update()

output = file.GetOutput()
scalar_range = output.GetScalarRange()

mapper = vtk.vtkDataSetMapper()
mapper.SetInputData(output)
mapper.SetScalarRange(scalar_range)

actor = vtk.vtkActor()
actor.SetMapper(mapper)

renderer = vtk.vtkRenderer()
renderer.AddActor(actor)
renderer.SetBackground(0, 0, 0)

renWindow = vtk.vtkRenderWindow()
renWindow.AddRenderer(renderer)

interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(renWindow)
interactor.Initialize()
interactor.Start()