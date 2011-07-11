import vtk


# Interpolate FEM onto imaging data structures
vtkExodusIIReader = vtk.vtkExodusIIReader()
vtkExodusIIReader.SetFileName("phantomMesh.e")
vtkExodusIIReader.ExodusModelMetadataOn()
vtkExodusIIReader.UpdateInformation()
vtkExodusIIReader.Update()
#print vtkExodusIIReader


GeneralTransform = vtk.vtkGeneralTransform()
GeneralTransform.RotateX(-90.0)
GeneralTransform.RotateY(0.0)
GeneralTransform.RotateZ(2.0)
GeneralTransform.Translate([0.039,0.031,0.0036])
GeneralTransform.Scale([1.,1.,1.])
GeneralTransform.Update()
print GeneralTransform.GetNumberOfConcatenatedTransforms()
matrix = GeneralTransform.GetConcatenatedTransform(0).GetMatrix()
print matrix

AffineTransform = vtk.vtkTransform()
AffineTransform.Translate([0.039,0.031,-0.001])
AffineTransform.RotateZ(2.0)
AffineTransform.RotateY(0.0)
AffineTransform.RotateX(-90.0)
AffineTransform.Scale([1.,1.,1.])
AffineTransform.Update()
print AffineTransform.GetNumberOfConcatenatedTransforms()
affinematrix = AffineTransform.GetConcatenatedTransform(0).GetMatrix()
print affinematrix
 
vec = [0.,0.,-1.]
vec1 = [0.,0.,.035]
print AffineTransform.TransformPoint(vec)
print AffineTransform.TransformVector(vec)
print AffineTransform.TransformPoint(vec1)
print AffineTransform.TransformVector(vec1)
dataImporter = vtk.vtkImageImport() 
#AppendBlocks = vtk.vtkAppendFilter()
#input = vtkExodusIIReader.GetOutput()
#idi = 0 
#if input.IsA("vtkMultiBlockDataSet"):
#    iter = input.NewIterator()
#    iter.UnRegister(None)
#    iter.InitTraversal()
#    while not iter.IsDoneWithTraversal():
#        curInput = iter.GetCurrentDataObject()
#        TransformFilter = vtk.vtkTransformFilter()
#        TransformFilter.SetTransform(GeneralTransform)
#        TransformFilter.SetInput( curInput )
#        TransformFilter.Update()
#        AppendBlocks.AddInput( TransformFilter.GetOutput() ) 
#        AppendBlocks.Update( ) 
#        #ExodusIIWrite = vtk.vtkExodusIIWriter()
#        #ExodusIIWrite.SetFileName("phantomMeshXform%02d.e" % idi)
#        #ExodusIIWrite.SetInput( TransformFilter.GetOutput() )
#        #ExodusIIWrite.SetModelMetadata( vtkExodusIIReader.GetExodusModel() )
#        #ExodusIIWrite.Update()
#        iter.GoToNextItem();
#        idi = idi + 1 
##
#vtkExodusIIWriter = vtk.vtkExodusIIWriter()
#vtkExodusIIWriter.SetFileName("phantomMeshXform.e")
#vtkExodusIIWriter.SetModelMetadata( vtkExodusIIReader.GetExodusModel() )
#vtkExodusIIWriter.SetInput( AppendBlocks.GetOutput() )
#vtkExodusIIWriter.Update()
