# Read DAKOTA parameters file (aprepro or standard format) and call a
# Read DAKOTA parameters file (aprepro or standard format) and call a
# Python module for fem analysis.

# DAKOTA will execute this script as
#   gPCWFS.py params.in results.out
# so sys.argv[1] will be the parameters file and
#    sys.argv[2] will be the results file to return to DAKOTA

# necessary python modules
import sys
import re
import os
import pickle
import ConfigParser

SolnOutputTemplate = "variable.%s.%04d"
MeshOutputTemplate = "fem_data.%04d.e" 

brainNekDIR = '/workarea/fuentes/braincode/tym1' 

setuprcFile = 'pigLiversetuprc'

workDirectory = 'optpp_pds'

##################################################################
def ComputeObjective(SEMDataDirectory,MRTIDirectory):
  import vtk
  import vtk.util.numpy_support as vtkNumPy 
  print "using vtk version", vtk.vtkVersion.GetVTKVersion()

  ObjectiveFunction = 0.0
  for timeID in [2]:
    # read SEM data
    vtkSEMReader = vtk.vtkXMLUnstructuredGridReader()
    vtufileName = "%s/%d.vtu" % (SEMDataDirectory,timeID)
    vtkSEMReader.SetFileName( vtufileName )
    vtkSEMReader.SetPointArrayStatus("Temperature",1)
    vtkSEMReader.Update()

    # load image 
    vtkImageReader = vtk.vtkDataSetReader() 
    vtkImageReader.SetFileName('%s/temperature.%04d.vtk' % (MRTIDirectory, timeID) )
    vtkImageReader.Update() 
    ## image_cells = vtkImageReader.GetOutput().GetPointData() 
    ## data_array = vtkNumPy.vtk_to_numpy(image_cells.GetArray('scalars')) 
    
    # extract voi for QOI
    vtkVOIExtract = vtk.vtkExtractVOI() 
    vtkVOIExtract.SetInput( vtkImageReader.GetOutput() ) 
    VOI = [10,100,100,150,0,0]
    vtkVOIExtract.SetVOI( VOI ) 
    vtkVOIExtract.Update()
    mrti_point_data= vtkVOIExtract.GetOutput().GetPointData() 
    mrti_array = vtkNumPy.vtk_to_numpy(mrti_point_data.GetArray('image_data')) 
    #print mrti_array
    #print type(mrti_array)

    # reuse ShiftScale Geometry
    vtkResample = vtk.vtkCompositeDataProbeFilter()
    vtkResample.SetInput( vtkImageReader.GetOutput() )
    vtkResample.SetSource( vtkVOIExtract.GetOutput() ) 
    vtkResample.Update()

    fem_point_data= vtkResample.GetOutput().GetPointData() 
    fem_array = vtkNumPy.vtk_to_numpy(fem_point_data.GetArray('Temperature')) 
    #print fem_array 
    #print type(fem_array )

    # accumulate objective function
    diff =  mrti_array-fem_array
    diffsq =  diff**2
    ObjectiveFunction = ObjectiveFunction + diffsq.sum()
  return ObjectiveFunction 
  
##################################################################
def brainNekWrapper(**kwargs):
  """
  call brainNek code 
  """
  outputMaterialFile = '%s/material_types.%04d.setup' % (workDirectory,kwargs['fileID'])
  print 'writing', outputMaterialFile 
  # open makefile
  fileHandle = file(outputMaterialFile   ,'w')
  fileHandle.write('[MATERIAL PROPERTIES]\n'  )
  fileHandle.write('# Name,      Type index, Density, Specific Heat, Conductivity, Perfusion, Absorption, Scattering, Anisotropy\n'  )

  print kwargs
  
  variableDictionary = kwargs['cv']
  fileHandle.write('Brain     0           1045     3640           %s        %s     %s      %s      %s \n' % ( variableDictionary['k_0_healthy'  ], variableDictionary['w_0_healthy'  ], variableDictionary['mu_s_healthy' ], variableDictionary['mu_a_healthy' ], variableDictionary['anfact'       ])
 )
  fileHandle.flush()
  fileHandle.close()

  # build command to run brainNek
  brainNekCommand = "%s/main ./%s -heattransfercoefficient %s -coolanttemperature  %s > %s/run.%04d.log" % (brainNekDIR , setuprcFile,variableDictionary['robin_coeff'  ], variableDictionary['probe_init'   ], workDirectory ,kwargs['fileID'])

  print brainNekCommand 
  os.system(brainNekCommand )
# end def brainNekWrapper:
##################################################################
def ParseInput(param_file):
  # ----------------------------
  # Parse DAKOTA parameters file
  # ----------------------------
  
  # setup regular expressions for parameter/label matching
  e = '-?(?:\\d+\\.?\\d*|\\.\\d+)[eEdD](?:\\+|-)?\\d+' # exponential notation
  f = '-?\\d+\\.\\d*|-?\\.\\d+'                        # floating point
  i = '-?\\d+'                                         # integer
  value = e+'|'+f+'|'+i                                # numeric field
  tag = '\\w+(?::\\w+)*'                               # text tag field
  
  # regular expression for aprepro parameters format
  aprepro_regex = re.compile('^\s*\{\s*(' + tag + ')\s*=\s*(' + value +')\s*\}$')
  # regular expression for standard parameters format
  standard_regex = re.compile('^\s*(' + value +')\s+(' + tag + ')$')
  
  # open DAKOTA parameters file for reading
  paramsfile = open(param_file, 'r')

  fileID = int(param_file.split(".").pop())
  #fileID = int(os.getcwd().split(".").pop())
  
  # extract the parameters from the file and store in a dictionary
  paramsdict = {}
  for line in paramsfile:
      m = aprepro_regex.match(line)
      if m:
          paramsdict[m.group(1)] = m.group(2)
      else:
          m = standard_regex.match(line)
          if m:
              paramsdict[m.group(2)] = m.group(1)
  
  paramsfile.close()
  
  # crude error checking; handle both standard and aprepro cases
  num_vars = 0
  if ('variables' in paramsdict):
      num_vars = int(paramsdict['variables'])
  elif ('DAKOTA_VARS' in paramsdict):
      num_vars = int(paramsdict['DAKOTA_VARS'])
  
  num_fns = 0
  if ('functions' in paramsdict):
      num_fns = int(paramsdict['functions'])
  elif ('DAKOTA_FNS' in paramsdict):
      num_fns = int(paramsdict['DAKOTA_FNS'])
  
  # initialize dictionary
  fem_params =  {} 

  # -------------------------------
  # Convert and send to application
  # -------------------------------
  
  # set up the data structures the rosenbrock analysis code expects
  # for this simple example, put all the variables into a single hardwired array
  continuous_vars = {} 

  try:
    continuous_vars['k_0_healthy' ] = paramsdict['k_0_healthy' ]
  except KeyError:
    pass

  try:
    continuous_vars['k_1'         ] = paramsdict['k_1'         ]
  except KeyError:
    pass

  try:
    continuous_vars['k_0_tumor'   ] = paramsdict['k_0_tumor'   ]
  except KeyError:
    pass

  try:
    continuous_vars['mu_a_healthy'] = paramsdict['mu_a_healthy']
  except KeyError:
    pass

  try:
    continuous_vars['mu_a_coag'  ]  = paramsdict['mu_a_coag'  ]
  except KeyError:
    pass

  try:
    continuous_vars['mu_a_tumor'  ] = paramsdict['mu_a_tumor'  ]
  except KeyError:
    pass

  try:
    continuous_vars['mu_s_healthy'] = paramsdict['mu_s_healthy']
  except KeyError:
    pass

  try:
    continuous_vars['mu_s_coag'  ] = paramsdict['mu_s_coag'  ]
  except KeyError:
    pass

  try:
    continuous_vars['mu_s_tumor'  ] = paramsdict['mu_s_tumor'  ]
  except KeyError:
    pass

  try:
    continuous_vars['w_0_healthy'] = paramsdict['w_0_healthy' ]  
  except KeyError:
    pass

  try:
    continuous_vars['w_1_coag'  ] = paramsdict['w_1_coag'   ] 
  except KeyError:
    pass
  
  try:
    continuous_vars['w_0_tumor'  ] = paramsdict['w_0_tumor'   ] 
  except KeyError:
    pass
  
  try:
    continuous_vars['anfact'] = paramsdict['anfact_healthy'   ] 
  except KeyError:
    pass
  
  try:
    continuous_vars['probe_init'] = paramsdict['probe_init'   ] 
  except KeyError:
    pass
  
  try:
    continuous_vars['robin_coeff'] = paramsdict['robin_coeff'   ] 
  except KeyError:
    pass
  
  try:
    active_set_vector = [ int(paramsdict['ASV_%d:response_fn_%d' % (i,i) ]) for i in range(1,num_fns+1)  ] 
  except KeyError:
    active_set_vector = [ int(paramsdict['ASV_%d:obj_fn' % (i) ]) for i in range(1,num_fns+1)  ] 
  
  # store dakota vars
  fem_params['cv']         = continuous_vars
  fem_params['asv']        = active_set_vector
  fem_params['functions']  = num_fns
  fem_params['fileID']     = fileID 

  return fem_params
  ## ----------------------------
  ## Return the results to DAKOTA
  ## ----------------------------
  #
  #if (fem_results['rank'] == 0 ):
  #  # write the results.out file for return to DAKOTA
  #  # this example only has a single function, so make some assumptions;
  #  # not processing DVV
  #  outfile = open('results.out.tmp.%d' % fileID, 'w')
  #  
  #  # write functions
  #  for func_ind in range(0, num_fns):
  #      if (active_set_vector[func_ind] & 1):
  #          functions = fem_results['fns']    
  #          outfile.write(str(functions[func_ind]) + ' f' + str(func_ind) + '\n')
  #  
  #  ## write gradients
  #  #for func_ind in range(0, num_fns):
  #  #    if (active_set_vector[func_ind] & 2):
  #  #        grad = rosen_results['fnGrads'][func_ind]
  #  #        outfile.write('[ ')
  #  #        for deriv in grad: 
  #  #            outfile.write(str(deriv) + ' ')
  #  #        outfile.write(']\n')
  #  #
  #  ## write Hessians
  #  #for func_ind in range(0, num_fns):
  #  #    if (active_set_vector[func_ind] & 4):
  #  #        hessian = rosen_results['fnHessians'][func_ind]
  #  #        outfile.write('[[ ')
  #  #        for hessrow in hessian:
  #  #            for hesscol in hessrow:
  #  #                outfile.write(str(hesscol) + ' ')
  #  #            outfile.write('\n')
  #  #        outfile.write(']]')
  #  #
  #  outfile.close();outfile.flush
  #  #
  #  ## move the temporary results file to the one DAKOTA expects
  #  #import shutil
  #  #shutil.move('results.out.tmp.%d' % fileID, sys.argv[2])
# end def ParseInput:
##################################################################

# setup command line parser to control execution
from optparse import OptionParser
parser = OptionParser()
parser.add_option( "--run_fem","--param_file", 
                  action="store", dest="param_file", default=None,
                  help="run code with parameter FILE", metavar="FILE")
(options, args) = parser.parse_args()


ComputeObjective("outputs/10-23-13-pigLiver-1","mrti")

if (options.param_file != None):
  fem_params = ParseInput(options.param_file)

  # FIXME link mesh directory
  linkDirectoryList = ['occa','libocca','meshes']
  for targetDirectory in linkDirectoryList:
    linkcommand = 'ln -sf %s/%s %s' % (brainNekDIR,targetDirectory ,targetDirectory  )
    print linkcommand 
    os.system(linkcommand )

  # execute the rosenbrock analysis as a separate Python module
  print "Running BrainNek..."
  fem_results = brainNekWrapper(**fem_params)
  print "BrainNek complete."

else:
  parser.print_help()
  print options
