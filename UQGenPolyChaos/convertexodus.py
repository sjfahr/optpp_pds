import os
import vtk
import sys

#datadirectory ="/workarea/fuentes/github/UQGenPolyChaos/biotex_090318_751642_treat_0"
datadirectory = sys.argv[1]
directoryList = os.listdir(datadirectory )
fileList      = filter(lambda x: x.split(".").pop() == "e" ,directoryList) 

for stats_file in fileList:       
  vtkExodusIIReader = vtk.vtkExodusIIReader()
  vtkExodusIIReader.SetFileName("%s/%s" % (datadirectory,stats_file))
  vtkExodusIIReader.Update()
  numberofresultarrays = vtkExodusIIReader.GetNumberOfPointResultArrays()
  print stats_file, numberofresultarrays
  for resultarrayindex in range(numberofresultarrays):
	resultarrayname = vtkExodusIIReader.GetPointResultArrayName(resultarrayindex)
	vtkExodusIIReader.SetPointResultArrayStatus( "%s" % (resultarrayname),1)
	#print resultarrayname
  vtkExodusIIReader.ExodusModelMetadataOn ()
  vtkExodusIIReader.Update()
  exodusObject = vtkExodusIIReader.GetOutput()
  ##print exodusObject.IsA("vtkMultiBlockDataSet")
  
  ## vtkXMLMultiBlockDataWriter = vtk.vtkXMLMultiBlockDataWriter()
  ## vtkXMLMultiBlockDataWriter.SetFileName("test.vtm")
  ## vtkXMLMultiBlockDataWriter.SetInput(exodusObject)
  ## vtkXMLMultiBlockDataWriter.Write()
  
  
  iter = exodusObject.NewIterator()
  iter.UnRegister(None)
  iter.InitTraversal()
  curInput = iter.GetCurrentDataObject()
  
  vtkUnstructuredGridWriter = vtk.vtkUnstructuredGridWriter()
  filebase = stats_file.split('.')
  filebase.pop()
  vtkUnstructuredGridWriter.SetFileName("%s/%s.vtk" % (datadirectory,'.'.join(filebase)))
  vtkUnstructuredGridWriter.SetInput(curInput)
  vtkUnstructuredGridWriter.Write()
