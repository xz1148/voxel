import vtk
import numpy as np

def ViewPoints(p, p_size):
    num_p = p.shape[0]
    points = vtk.vtkPoints()
    vertices = vtk.vtkCellArray()
    vertices.InsertNextCell(num_p)
    for x in p:
        id = points.InsertNextPoint(x)
        vertices.InsertCellPoint(id)



    point = vtk.vtkPolyData()
    point.SetPoints(points)
    point.SetVerts(vertices)
    mapper = vtk.vtkPolyDataMapper()
    if vtk.VTK_MAJOR_VERSION <=5:
        mapper.SetInput(point)
    else:
        mapper.SetInputData(point)


    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetPointSize(p_size)



    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    actor.GetProperty().SetColor(1,0.5,1)
    ren.AddActor(actor)

    iren.Initialize()
    renWin.Render()
    iren.Start()

