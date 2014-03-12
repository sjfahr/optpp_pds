# necessary python modules
import sys
import re
import os
import scipy
import numpy

# FIXME global vars are prob bad idea...
ntime = 60
nsubstep = 1
acquisitionTime = 5.00
deltat = acquisitionTime / nsubstep
#FEMMeshFileName="/mnt/data/scratch/fuentes/data/mdacc/uq_canine_feb07/meshTemplateLowerRes.e"
#FEMMeshFileName="/work/00131/fuentes/data/mdacc/uq_canine_feb07/meshTemplateQuarter2.e"
#FEMMeshFileName="/data/fuentes/mdacc/uq_canine_feb07/meshTemplateQuarter1.e"
#FEMMeshFileName="/data/cjmaclellan/mdacc/nano/phantomMeshFullTess.e"
FEMMeshFileName="/data/cjmaclellan/mdacc/deltap_phantom_oct10/phantomMeshFull.e"
SolnOutputTemplate = "soln.%04d.out"
#print FEMMeshFileName

def ProjectImagingMesh(fem_mesh_file): 
  import vtk
  import vtk.util.numpy_support as vtkNumPy
  import scipy.io as scipyio
  # echo vtk version info
  print "using vtk version", vtk.vtkVersion.GetVTKVersion()
  # FIXME notice that order of operations is IMPORTANT
  # FIXME translation followed by rotation will give different results
  # FIXME than rotation followed by translation
  # FIXME Translate -> RotateZ -> RotateY -> RotateX -> Scale seems to be the order of paraview
  AffineTransform = vtk.vtkTransform()
  # should be in meters
  AffineTransform.Translate(.0000001,.0000001,.0000001)
  AffineTransform.RotateZ( 0 )
  AffineTransform.RotateY( 0 )
  AffineTransform.RotateX( 0 )
  AffineTransform.Scale([1.,1.,1.])
  # get homogenius 4x4 matrix of the form
  # A | b
  # matrix = -----
  # 0 | 1
  #
  matrix = AffineTransform.GetConcatenatedTransform(0).GetMatrix()
  #print matrix
  RotationMatrix = [[matrix.GetElement(0,0),matrix.GetElement(0,1),matrix.GetElement(0,2)],
                    [matrix.GetElement(1,0),matrix.GetElement(1,1),matrix.GetElement(1,2)],
                    [matrix.GetElement(2,0),matrix.GetElement(2,1),matrix.GetElement(2,2)]]
  Translation = [matrix.GetElement(0,3),matrix.GetElement(1,3),matrix.GetElement(2,3)]
  #print RotationMatrix ,Translation
  

  #laserTip = AffineTransform.TransformPoint( laserTip )
  #laserOrientation = AffineTransform.TransformVector( laserOrientation )

  # read imaging data geometry that will be used to project FEM data onto
  #vtkReader = vtk.vtkXMLImageDataReader()
  vtkReader = vtk.vtkDataSetReader()
  #vtkReader.SetFileName('/data/cjmaclellan/mdacc/nano/spio/spioVTK/67_11/tmap_67_11.0000.vtk' )
  vtkReader.SetFileName('/data/cjmaclellan/mdacc/deltap_phantom_oct10/VTKtmaps/S695/S695.0000.vtk' )
  #vtkReader.SetFileName('/data/cjmaclellan/mdacc/nano/nrtmapsVTK/R695/R695.0000.vtk' ) 
  vtkReader.Update()
  templateImage = vtkReader.GetOutput()
  dimensions = templateImage.GetDimensions()
  spacing = templateImage.GetSpacing()
  origin = templateImage.GetOrigin()
  print spacing, origin, dimensions 
  #fem.SetImagingDimensions( dimensions ,origin,spacing)

  ## project imaging onto fem mesh
  #if (vtk != None):
  # vtkImageReader = vtk.vtkDataSetReader()
  # vtkImageReader.SetFileName('/data/fuentes/mdacc/uqModelStudy/dog1_980/temperature.%04d.vtk' % timeID )
  # vtkImageReader.Update()
  # image_cells = vtkImageReader.GetOutput().GetPointData()
  # data_array = vtkNumPy.vtk_to_numpy(image_cells.GetArray('scalars'))
  # v1 = PETSc.Vec().createWithArray(data_array, comm=PETSc.COMM_SELF)
  # fem.ProjectImagingToFEMMesh("MRTI", v1)
  # fem.StoreSystemTimeStep("MRTI",timeID )
  #
  # # extract voi for QOI
  # vtkVOIExtract = vtk.vtkExtractVOI()
  # vtkVOIExtract.SetInput( vtkImageReader.GetOutput() )
  # VOI = [10,100,100,150,0,0]
  # vtkVOIExtract.SetVOI( VOI )
  # vtkVOIExtract.Update()
  # mrti_point_data= vtkVOIExtract.GetOutput().GetPointData()
  # mrti_array = vtkNumPy.vtk_to_numpy(mrti_point_data.GetArray('scalars'))
  # #print mrti_array
  # #print type(mrti_array)

  # Interpolate FEM onto imaging data structures
  vtkExodusIIReader = vtk.vtkExodusIIReader()
  vtkExodusIIReader.SetFileName( fem_mesh_file )
  vtkExodusIIReader.SetPointResultArrayStatus("u0",1)
  vtkExodusIIReader.SetPointResultArrayStatus("u0*",1)
  vtkExodusIIReader.SetPointResultArrayStatus("u1",1)

  matsize = int(dimensions[1])#get matrix size
  #preallocate size of arrays
  u0_array_1 = scipy.zeros((matsize,matsize))
  u0_array_2 = scipy.zeros((matsize,matsize,ntime*nsubstep))
  u1_array_1 = scipy.zeros((matsize,matsize))
  u1_array_2 = scipy.zeros((matsize,matsize,ntime*nsubstep))
  u0star_array_1 = scipy.zeros((matsize,matsize))
  u0star_array_2 = scipy.zeros((matsize,matsize,ntime*nsubstep))

  #for timeID in range(1,2):
  for timeID in range(1,ntime*nsubstep):
    vtkExodusIIReader.SetTimeStep(timeID-1)
    vtkExodusIIReader.Update()
    
    # apply the transform
    TransformedFEMMesh = None
    if vtkExodusIIReader.GetOutput().IsA("vtkMultiBlockDataSet"):
      AppendBlocks = vtk.vtkAppendFilter()
      iter = vtkExodusIIReader.GetOutput().NewIterator()
      iter.UnRegister(None)
      iter.InitTraversal()
      # loop over blocks...
      while not iter.IsDoneWithTraversal():
          curInput = iter.GetCurrentDataObject()
          vtkTransformFEMMesh = vtk.vtkTransformFilter()
          vtkTransformFEMMesh.SetTransform( AffineTransform )
          vtkTransformFEMMesh.SetInput( curInput )
          vtkTransformFEMMesh.Update()
          AppendBlocks.AddInput( vtkTransformFEMMesh.GetOutput() )
          AppendBlocks.Update( )
          iter.GoToNextItem();
      TransformedFEMMesh = AppendBlocks.GetOutput()
    else:
      vtkTransformFEMMesh = vtk.vtkTransformFilter()
      vtkTransformFEMMesh.SetTransform( AffineTransform )
      vtkTransformFEMMesh.SetInput( vtkExodusIIReader.GetOutput() )
      vtkTransformFEMMesh.Update()
      TransformedFEMMesh = vtkTransformFEMMesh.GetOutput()

    # reflect
    #vtkReflectX = vtk.vtkReflectionFilter()
    #vtkReflectX.SetPlaneToXMin()
    #vtkReflectX.SetInput( TransformedFEMMesh )
    #vtkReflectX.Update()

    # reflect
    #vtkReflectZ = vtk.vtkReflectionFilter()
    #vtkReflectZ.SetPlaneToZMax()
    #vtkReflectZ.SetInput( vtkReflectX.GetOutput() )
    #vtkReflectZ.Update()

    # reuse ShiftScale Geometry
    vtkResample = vtk.vtkCompositeDataProbeFilter()
    vtkResample.SetInput( templateImage )
    vtkResample.SetSource( TransformedFEMMesh )
    #vtkResample.SetSource( vtkReflectZ.GetOutput() )
    vtkResample.Update()
    fem_point_data= vtkResample.GetOutput().GetPointData()
    u0_array = vtkNumPy.vtk_to_numpy(fem_point_data.GetArray('u0'))
    u0star_array = vtkNumPy.vtk_to_numpy(fem_point_data.GetArray('u0*'))
    u1_array = vtkNumPy.vtk_to_numpy(fem_point_data.GetArray('u1'))
 
    
    #go from 1x256^2 array to 256X256 array for each timestep
    for nn in range(0,matsize):
     u0_array_1[nn,:]=u0_array[nn*matsize:(nn+1)*matsize]
     u0star_array_1[nn,:]=u0star_array[nn*matsize:(nn+1)*matsize]
     u1_array_1[nn,:]=u1_array[nn*matsize:(nn+1)*matsize]
   
    #apply proper rotations/reflections and combine 2D arrays into 3D array 
    u0_array_2[:,:,timeID-1]=numpy.fliplr(numpy.rot90(u0_array_1,k=3))
    u0star_array_2[:,:,timeID-1]=numpy.fliplr(numpy.rot90(u0star_array_1,k=3))
    u1_array_2[:,:,timeID-1]=numpy.fliplr(numpy.rot90(u1_array_1,k=3))
    
    # write numpy to disk in matlab
    #scipyio.savemat("MS795.%04d.mat" % (timeID), {'u0':u1_array,'MRTI0':MRTI0_array })

    # write output
    print "writing ", timeID
    vtkStatsWriter = vtk.vtkDataSetWriter()
    vtkStatsWriter.SetFileTypeToBinary()
    vtkStatsWriter.SetFileName("test.%04d.vtk" % timeID )
    vtkStatsWriter.SetInput(vtkResample.GetOutput())
    vtkStatsWriter.Update()

  scipyio.savemat("S695.mat",{'ModelFluence':u1_array_2,'MRTI':u0star_array_2,'ModelTemp':u0_array_2})
  
#ProjectImagingMesh('/data/cjmaclellan/mdacc/nano/spio/spioVTK/335_11/fem_data_335_11.e')
#ProjectImagingMesh('/data/cjmaclellan/mdacc/nano/nrtmapsVTK/R695/fem_data_NR.e')
ProjectImagingMesh('/data/cjmaclellan/mdacc/deltap_phantom_oct10/fem_data.0001.e')
#ProjectImagingMesh(FEMMeshFileName)
