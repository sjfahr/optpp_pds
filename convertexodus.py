import os
import vtk

datadirectory ="/workarea/fuentes/github/UQGenPolyChaos/biotex_090318_751642_treat_0"
directoryList = os.listdir(datadirectory )
fileList      = filter(lambda x: x.split(".").pop(0) == "fem_stats" ,directoryList) 

for stats_file in fileList:       
  vtkExodusIIReader = vtk.vtkExodusIIReader()
  vtkExodusIIReader.SetFileName("%s/%s" % (datadirectory,stats_file))
  vtkExodusIIReader.SetPointResultArrayStatus("Varu0Mean",1)
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
