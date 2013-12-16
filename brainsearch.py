# Read DAKOTA parameters file (aprepro or standard format) and call a
# Python module for fem analysis.
# DAKOTA will execute this script 

# necessary python modules
import sys
import re
import os
import scipy.io as scipyio
import ConfigParser

brainNekDIR     = '/workarea/fuentes/braincode/tym1' 
workDirectory   = 'optpp_pds'
outputDirectory = '/dev/shm/outputs/dakota/%04d'
outputDirectory = '/tmp/outputs/dakota/%04d'

# database and run directory have the same structure
databaseDIR     = 'database/'
# $ ls database workdir/
# database:
# Patient0002/  Patient0003/  Patient0004/  Patient0005/  Patient0006/  Patient0007/  Patient0008/
# 
# workdir/:
# Patient0002/  Patient0003/  Patient0004/  Patient0005/  Patient0006/  Patient0007/  Patient0008/
# $ ls database/Patient0002/ workdir/Patient0002
# database/Patient0002/:
# 000/  001/  010/  017/  020/  021/  022/
# 
# workdir/Patient0002:
# 000/  001/  010/  017/  020/  021/  022/

caseFunctionTemplate = \
"""
// Prototypes
occaDeviceFunction datafloat laserPower(datafloat time);
occaDeviceFunction datafloat initialTemperature(datafloat x, datafloat y, datafloat z);
occaDeviceFunction datafloat sourceFunction(datafloat x , datafloat y , datafloat z, 
			 datafloat x0, datafloat y0, datafloat z0, 
			 datafloat volumeFraction, 
			 datafloat muA_muTr      , datafloat muEff);
occaDeviceFunction datafloat DirichletTemp(unsigned int bcTag, datafloat time,
			datafloat x, datafloat y, datafloat z);
occaDeviceFunction datafloat RobinCoeff(unsigned int bcTag, datafloat x, datafloat y, datafloat z,
		     datafloat kappa, datafloat h);
occaDeviceFunction datafloat NeumannDeriv(unsigned int bcTag, datafloat time, 
		       datafloat x, datafloat y, datafloat z,
		       datafloat nx, datafloat ny, datafloat nz);
occaDeviceFunction datafloat exactSolution(datafloat x, datafloat y, datafloat z, datafloat time,
			datafloat kappa, datafloat lambda);

/*
 * Compile-time definitions
 *   - bodyTemperature    = ambient body temperature
 *   - coolantTemperature = probe coolant temperature
 *   - laserMaxPower      = reference laser power
 */

/// Laser power as a function of time
/**
 * @param time
 */
occaDeviceFunction datafloat laserPower(datafloat time) {
%s
}

/// Initial temperature
/**
 * Boundary conditions will be enforced afterwards
 * @param x
 * @param y
 * @param z
 * @param bodyTemperature
 * @return initial temperature
 */
occaDeviceFunction datafloat initialTemperature(datafloat x, datafloat y, datafloat z) {
  return bodyTemperature;
}

/// Heating at a point due to a region of the laser tip
/**
 * @param x
 * @param y
 * @param z
 * @param x0 x-coordinate of centroid of laser tip region
 * @param y0 y-coordinate of centroid of laser tip region
 * @param z0 z-coordinate of centroid of laser tip region
 * @param volumeFraction volume fraction of laser tip region relative to the 
 *          entire laser tip
 * @param mu_a absorption coefficient of laser light in tissue
 * @param mu_eff effective absorption (\f$\mu_\text{eff}=\sqrt{3\mu_a\mu_{tr}}\f$)
 * @param mu_tr transport coefficient (\f$\mu_{tr}=\mu_a + \mu_s (1-g)\f$)
 * @return contribution of source point to heating function
 */
occaDeviceFunction datafloat sourceFunction(datafloat x , datafloat y , datafloat z, 
			 datafloat x0, datafloat y0, datafloat z0, 
			 datafloat volumeFraction, 
			 datafloat muA_muTr      , datafloat muEff) {
  // Distance between point and source point
  datafloat dist = (x - x0)*(x - x0) + (y - y0)*(y - y0) + (z - z0)*(z - z0);
  dist = sqrt(dist);

  // Choose minimum distance to avoid dividing by zero
  if(dist < 1e-6)
    return 0;

  // Return contribution to forcing function
  return 0.75*M_1_PI*muA_muTr*volumeFraction*exp(-muEff*dist)/dist;
}


/// Returns the temperature corresponding to a Dirichlet boundary condition
/**
 * @param bcTag type of boundary condition
 *          - 1 = body temperature Dirichlet boundary condition
 *          - 2 = coolant temperature Dirichlet boundary condition
 * @param x
 * @param y
 * @param z
 * @param time
 * @return Dirichlet boundary condition temperature
 */
occaDeviceFunction datafloat DirichletTemp(unsigned int bcTag, datafloat time, 
			datafloat x, datafloat y, datafloat z) {
  switch(bcTag) {
  case 1:  return bodyTemperature;
  case 2:  return coolantTemperature;
  default: break;
  }
  
  return bodyTemperature;
}

/// Returns the coefficient corresponding to a Robin boundary condition
/**
 * We assume a Robin boundary condition of the form 
 *   \f[\kappa\frac{\partial u}{\partial n}=-\alpha\left(u-u_b\right)\f]
 * @param bcTag type of boundary condition
 *          - 3 = Neumann boundary condition (\f$\alpha=0\f$)
 *          - 4 = Robin condition at probe
 * @param x
 * @param y
 * @param z
 * @param kappa thermal conductivity
 * @param h heat transfer coefficient
 * @return \f$\alpha\f$
 */
occaDeviceFunction datafloat RobinCoeff(unsigned int bcTag, 
		     datafloat x, datafloat y, datafloat z,
		     datafloat kappa, datafloat h) {
  switch(bcTag) {
  case 3:
    return 0;
  case 4:
    return h; // Heat transfer coefficient
  default: return 0;
  }
}		  


/// Returns the derivative corresponding to a Neumann boundary condition
/**
 * Note: not currently used
 * @param bcTag type of boundary condition
 * @param time
 * @param x
 * @param y
 * @param z
 * @param nx x-coordinate of surface normal vector
 * @param ny y-coordinate of surface normal vector
 * @param nz z-coordinate of surface normal vector
 */
occaDeviceFunction datafloat NeumannDeriv(unsigned int bcTag, datafloat time, 
		       datafloat x, datafloat y, datafloat z,
		       datafloat nx, datafloat ny, datafloat nz){
  // Homogeneous Neumann
  return 0;
}

/// Analytic solution
/**
 * Note: an analytic solution is not known for this case. 
 * The numerical solution can be compared with the analytic solution, if it 
 *   is known
 * @param kappa tissue thermal conductivity
 * @param lambda (tissue density)*(tissue specific heat)/dt 
 *                 + (perfusion)*(blood specific heat)
 * @param x
 * @param y
 * @param z
 * @param time
 * @return temperature
 */
occaDeviceFunction datafloat exactSolution(datafloat x, datafloat y, datafloat z, datafloat time,
			datafloat kappa, datafloat lambda){
  return bodyTemperature;
}
"""

setuprcTemplate = \
"""
[THREAD MODEL]
OpenCL

[CASE FILE]
%s/case.%04d.setup

[MESH FILE]
meshes/cooledConformMesh.inp

[MRI FILE]
./mridata.setup

[POLYNOMIAL ORDER]
3

[DT]
0.25

[FINAL TIME]
%f

[PCG TOLERANCE]
1e-6

[GPU PLATFORM]
0

[GPU DEVICE]
1

[SCREENSHOT OUTPUT]
%s

[SCREENSHOT INTERVAL]
%f
"""

caseFileTemplate = \
"""
# case setup file
# all physical quantities should be in MKS units and degrees Celsius

# Data is for porcine liver (Roggan and Muller 1994)

[FUNCTION FILE]
%s/casefunctions.%04d.occa

[HAS EXACT SOLUTION]
0

# first row is center of probe, second row specifies direction
[PROBE POSITION]
0 0 0
0 0 1

[LASER MAXIMUM POWER]
15

[BODY TEMPERATURE]
%s

[BLOOD TEMPERATURE]
%s

[COOLANT TEMPERATURE]
%s

[BLOOD SPECIFIC HEAT]
3840

[DAMAGE FREQUENCY FACTOR]
1e70

[DAMAGE ACTIVATION ENERGY]
4e5

[GAS CONSTANT]
8.314

[PROBE HEAT TRANSFER COEFFICIENT]
0

[MESH BLOCKS - COMPUTATIONAL DOMAIN]
brain

[MESH BLOCKS - LASER TIP]
laserTip

[MESH BLOCKS - DIRICHLET (BODY TEMPERATURE)]

[MESH BLOCKS - DIRICHLET (COOLANT TEMPERATURE)]

[MESH SIDE SETS - DIRICHLET (BODY TEMPERATURE)]
regionBoundary

[MESH SIDE SETS - DIRICHLET (COOLANT TEMPERATURE)]

[MESH SIDE SETS - NEUMANN]
probeSurface

[MESH SIDE SETS - ROBIN]

[BRAIN TYPES - DIRICHLET (BODY TEMPERATURE)]

[BRAIN MATERIAL DATA FILE]
./material_data.setup

[BRAIN MATERIAL PROPERTIES FILE]
%s/material_types.%04d.setup

# Currently has material properties of water
[PROBE MATERIAL PROPERTIES]
# Name,   Density, Specific Heat, Conductivity, Absorption, Scattering, Anisotropy
catheter  1.0      4180           0.5985        500         14000       0.88
laserTip  1.0      4180           0.5985        500         14000       0.88
"""

##################################################################
def ComputeObjective(**kwargs):
  ObjectiveFunction = 0.0
  # initialize brainNek
  import brainNekLibrary
  import numpy
  # setuprc file
  outputSetupRCFile = '%s/setuprc.%04d' % (workDirectory,kwargs['fileID'])
  setup = brainNekLibrary.PySetupAide(outputSetupRCFile )
  brainNek = brainNekLibrary.PyBrain3d(setup);

  # write out initial mesh
  
  # FIXME vtk needs to be loaded AFTER kernel is built
  import vtk
  import vtk.util.numpy_support as vtkNumPy 
  print "using vtk version", vtk.vtkVersion.GetVTKVersion()

  # setup vtkUnstructuredGrid
  hexahedronGrid   = vtk.vtkUnstructuredGrid()
  numPoints = brainNek.GetNumberOfNodes( ) 
  numElems  = brainNek.GetNumberOfElements( ) 
  # initialize nodes and connectivity
  bNekNodes         = numpy.zeros(numPoints * 3,dtype=numpy.float32)
  bNekConnectivity  = numpy.zeros(numElems  * 8,dtype=numpy.int32)
  print "setting up hex mesh with %d nodes %d elem"  % (numPoints,numElems)

  # get nodes and connectivity from brainnek
  brainNek.GetNodes(   bNekNodes)       ;
  brainNek.GetElements(bNekConnectivity);
  # reshape for convience
  bNekNodes        = bNekNodes.reshape(      numPoints , 3)
  bNekConnectivity = bNekConnectivity.reshape(numElems , 8)

  # setup points
  hexahedronPoints = vtk.vtkPoints()
  for nodecoord in bNekNodes:
    ipoint = hexahedronPoints.InsertNextPoint(nodecoord );
  hexahedronGrid.SetPoints(hexahedronPoints);

  # setup elements
  for ielem in range(numElems):
    aHexahedron = vtk.vtkHexahedron()
    aHexahedron.GetPointIds().SetId(0,bNekConnectivity[ielem][0])
    aHexahedron.GetPointIds().SetId(1,bNekConnectivity[ielem][1])
    aHexahedron.GetPointIds().SetId(2,bNekConnectivity[ielem][2])
    aHexahedron.GetPointIds().SetId(3,bNekConnectivity[ielem][3])
    aHexahedron.GetPointIds().SetId(4,bNekConnectivity[ielem][4])
    aHexahedron.GetPointIds().SetId(5,bNekConnectivity[ielem][5])
    aHexahedron.GetPointIds().SetId(6,bNekConnectivity[ielem][6])
    aHexahedron.GetPointIds().SetId(7,bNekConnectivity[ielem][7])
    hexahedronGrid.InsertNextCell(aHexahedron.GetCellType(),
                                  aHexahedron.GetPointIds())
  print "done setting hex mesh with %d nodes %d elem"  % (numPoints,numElems)

  # setup solution
  bNekSoln = numpy.zeros(numPoints,dtype=numpy.float32)
  brainNek.getHostTemperature(bNekSoln )
  DeepCopy = 1
  vtkScalarArray = vtkNumPy.numpy_to_vtk( bNekSoln, DeepCopy) 
  vtkScalarArray.SetName("bioheat") 
  hexahedronGrid.GetPointData().SetScalars(vtkScalarArray);

  ## # dbg 
  ## brainNek.screenshot( 0.0 )

  # get registration parameters
  variableDictionary = kwargs['cv']

  # register the SEM data to MRTI
  AffineTransform = vtk.vtkTransform()
  AffineTransform.Translate([ 
    float(variableDictionary['x_displace']),
    float(variableDictionary['y_displace']),
    float(variableDictionary['z_displace'])
                            ])
  # FIXME  notice that order of operations is IMPORTANT
  # FIXME   translation followed by rotation will give different results
  # FIXME   than rotation followed by translation
  # FIXME  Translate -> RotateZ -> RotateY -> RotateX -> Scale seems to be the order of paraview
  AffineTransform.RotateZ( float(variableDictionary['z_rotate'  ] ) ) 
  AffineTransform.RotateY( float(variableDictionary['y_rotate'  ] ) )
  AffineTransform.RotateX( float(variableDictionary['x_rotate'  ] ) )
  AffineTransform.Scale([1.e0,1.e0,1.e0])
  SEMRegister = vtk.vtkTransformFilter()
  SEMRegister.SetInput( hexahedronGrid )
  SEMRegister.SetTransform(AffineTransform)
  SEMRegister.Update()

  ## vtkSEMReader = vtk.vtkXMLUnstructuredGridReader()
  ## SEMDataDirectory = outputDirectory % kwargs['UID']
  ## SEMtimeID = 0 
  ## vtufileName = "%s/%d.vtu" % (SEMDataDirectory,SEMtimeID)
  ## print "reading ", vtufileName 
  ## vtkSEMReader.SetFileName( vtufileName )
  ## vtkSEMReader.SetPointArrayStatus("Temperature",1)
  ## vtkSEMReader.Update()
  ## fem_point_data= vtkSEMReader.GetOutput().GetPointData() 
  ## tmparray = vtkNumPy.vtk_to_numpy(fem_point_data.GetArray('Temperature')) 

  # loop over time steps
  tstep = 0
  currentTime = 0.0

  # setup MRTI data read
  MRTItimeID  = 0
  MRTIInterval = 5.0

  # setup screen shot interval 
  screenshotNum = 1;
  screenshotTol = 1e-10;
  screenshotInterval = MRTIInterval ;

  ## # loop over time points of interest
  ## for (SEMtimeID,MRTItimeID) in [(1,kwargs['maxheat'])]:

  while( brainNek.timeStep(tstep * .25 ) ) :
    tstep = tstep + 1
    currentTime = tstep * .25
    ## if(currentTime+screenshotTol >= screenshotNum*screenshotInterval):
    ##    brainNek.getHostTemperature( bNekSoln )
    ##    screenshotNum = screenshotNum + 1;
    ##    print "get host data",bNekSoln 

    if(currentTime+screenshotTol >= MRTItimeID * MRTIInterval):
      # load image 
      mrtifilename = '%s/temperature.%04d.vtk' % (kwargs['mrti'], MRTItimeID) 
      print 'opening' , mrtifilename 
      vtkImageReader = vtk.vtkDataSetReader() 
      vtkImageReader.SetFileName(mrtifilename )
      vtkImageReader.Update() 
      ## image_cells = vtkImageReader.GetOutput().GetPointData() 
      ## data_array = vtkNumPy.vtk_to_numpy(image_cells.GetArray('scalars')) 
      
      # extract voi for QOI
      vtkVOIExtract = vtk.vtkExtractVOI() 
      vtkVOIExtract.SetInput( vtkImageReader.GetOutput() ) 
      vtkVOIExtract.SetVOI( kwargs['voi'] ) 
      vtkVOIExtract.Update()
      mrti_point_data= vtkVOIExtract.GetOutput().GetPointData() 
      mrti_array = vtkNumPy.vtk_to_numpy(mrti_point_data.GetArray('image_data')) 
      #print mrti_array
      #print type(mrti_array)

      # get brainNek solution 
      brainNek.getHostTemperature( bNekSoln )
      vtkScalarArray = vtkNumPy.numpy_to_vtk( bNekSoln, DeepCopy) 
      vtkScalarArray.SetName("bioheat") 
      hexahedronGrid.GetPointData().SetScalars(vtkScalarArray);

      # project SEM onto MRTI for comparison
      print 'resampling' 
      vtkResample = vtk.vtkCompositeDataProbeFilter()
      vtkResample.SetSource( SEMRegister.GetOutput() )
      vtkResample.SetInput( vtkVOIExtract.GetOutput() ) 
      vtkResample.Update()

      fem_point_data= vtkResample.GetOutput().GetPointData() 
      fem_array = vtkNumPy.vtk_to_numpy(fem_point_data.GetArray('bioheat')) 
      print 'resampled' ,fem_array 
      #print fem_array 
      #print type(fem_array )

      # FIXME  should this be different ?  
      SEMDataDirectory = outputDirectory % kwargs['UID']

      DebugObjective = True
      DebugObjective = False
      # write output
      if ( DebugObjective ):
        vtkSEMWriter = vtk.vtkXMLUnstructuredGridWriter()
        semfileName = "%s/semtransform%04d.vtu" % (SEMDataDirectory,MRTItimeID)
        print "writing ", semfileName 
        vtkSEMWriter.SetFileName( semfileName )
        vtkSEMWriter.SetInput(SEMRegister.GetOutput())
        #vtkSEMWriter.SetDataModeToAscii()
        vtkSEMWriter.Update()

      # write output
      if ( DebugObjective ):
         vtkTemperatureWriter = vtk.vtkDataSetWriter()
         vtkTemperatureWriter.SetFileTypeToBinary()
         roifileName = "%s/roi%04d.vtk" % (SEMDataDirectory,MRTItimeID)
         print "writing ", roifileName 
         vtkTemperatureWriter.SetFileName( roifileName )
         vtkTemperatureWriter.SetInput(vtkResample.GetOutput())
         vtkTemperatureWriter.Update()

      # accumulate objective function
      diff =  numpy.abs(mrti_array-fem_array)
      diffsq =  diff**2
      ObjectiveFunction = ObjectiveFunction + diff.sum()

      # update counter
      MRTItimeID = MRTItimeID + 1;

  return ObjectiveFunction 
# end def ComputeObjective:
##################################################################
def brainNekWrapper(**kwargs):
  """
  call brainNek code 
  """
  import math
  # occa case file
  outputOccaCaseFile = '%s/casefunctions.%04d.occa' % (workDirectory,kwargs['fileID'])
  print 'writing', outputOccaCaseFile 
  with file(outputOccaCaseFile, 'w') as occaCaseFileName: occaCaseFileName.write(caseFunctionTemplate % kwargs['ccode'] )

  # setuprc file
  outputSetupRCFile = '%s/setuprc.%04d' % (workDirectory,kwargs['fileID'])
  print 'writing', outputSetupRCFile 
  fileHandle = file(outputSetupRCFile ,'w')
  semwritetime = kwargs['semwritetime']
  # make sure write directory exists
  os.system('mkdir -p %s' % outputDirectory % kwargs['UID'] )
  fileHandle.write(setuprcTemplate % (workDirectory,kwargs['fileID'] ,semwritetime, outputDirectory % kwargs['UID'] ,semwritetime  ) )
  fileHandle.flush(); fileHandle.close()

  # get variables
  variableDictionary = kwargs['cv']

  # case file
  outputCaseFile = '%s/case.%04d.setup' % (workDirectory,kwargs['fileID'])
  print 'writing', outputCaseFile 
  with file(outputCaseFile , 'w') as fileHandle: fileHandle.write(caseFileTemplate % (workDirectory,kwargs['fileID'],variableDictionary['body_temp'],variableDictionary['body_temp'],variableDictionary['probe_init'],workDirectory,kwargs['fileID'])  )

#      mu_a_min               <      mu_a + (1-g) mu_s < mu_a_max + (1-g_min) mu_s_max
#         5.e-1               <          mu_tr         < 600. + .3 * 50000. 
#
#  sqrt( 3 * 5.e-1 * 5.e-1 )  <  sqrt( 3 mu_a  mu_tr ) < sqrt( 3 * 600. * (600. + .3 * 50000.) ) 
#  sqrt( 3 * 5.e-1 * 5.e-1 )  <        mu_eff          < sqrt( 3 * 600. * (600. + .3 * 50000.) ) 
#            8.e-1            <        mu_eff          <    5.3e3

  # materials
  outputMaterialFile = '%s/material_types.%04d.setup' % (workDirectory,kwargs['fileID'])
  print 'writing', outputMaterialFile 
  fileHandle = file(outputMaterialFile   ,'w')
  fileHandle.write('[MATERIAL PROPERTIES]\n'  )
  fileHandle.write('# Name,      Type index, Density, Specific Heat, Conductivity, Perfusion, Absorption, Scattering, Anisotropy\n'  )
  mu_s   = 8.e3
  anfact = .9
  mu_s_p = mu_s * (1.-anfact) 
  mu_eff = float(variableDictionary['mu_eff_healthy'])
  mu_a   =  0.5*( -mu_s_p + math.sqrt( mu_s_p * mu_s_p  + 4. * mu_eff * mu_eff  /3. ) )
  fileHandle.write('Brain     0           1045     3640           %s        %s     %f      %f      %f \n' % ( variableDictionary['k_0_healthy'  ], variableDictionary['w_0_healthy'  ], mu_a, mu_s, anfact )
 )
  fileHandle.flush(); fileHandle.close()

  ## # build command to run brainNek
  ## brainNekCommand = "%s/main %s -heattransfercoefficient %s -coolanttemperature  %s > %s/run.%04d.log 2>&1 " % (brainNekDIR , outputSetupRCFile ,variableDictionary['robin_coeff'  ], variableDictionary['probe_init'   ], workDirectory ,kwargs['fileID'])

  ## # system call to run brain code
  ## print brainNekCommand 
  ## os.system(brainNekCommand )
# end def brainNekWrapper:
##################################################################
def ParseInput(paramfilename):
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
  paramsfile = open(paramfilename, 'r')

  fileID = int(paramfilename.split(".").pop())
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

  DescriptorList = ['robin_coeff','probe_init','mu_eff_healthy','body_temp','anfact_healthy', 'mu_a_healthy','mu_s_healthy','k_0_healthy','w_0_healthy','x_displace','y_displace','z_displace','x_rotate','y_rotate','z_rotate']
  for paramname in DescriptorList:
    try:
      continuous_vars[paramname  ] = paramsdict[paramname ]
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
  fem_params['UID']        = int(paramfilename.split('/').pop(3))

  # parse file path
  locatemrti = paramfilename.split('/')
  locatemrti.pop()

  # database and run directory have the same structure
  fem_params['mrti']       = '%s/%s/%s/vtk/referenceBased/' % (databaseDIR,locatemrti[2],locatemrti[3])

  # get power file name
  inisetupfile  = "/".join(locatemrti)+"/setup.ini"
  config = ConfigParser.SafeConfigParser({})
  config.read(inisetupfile)
  fem_params['ccode']        = config.get('power','ccode')
  fem_params['semwritetime'] = config.getfloat('mrti','deltat') * config.getfloat('mrti','maxheat')
  fem_params['maxheat']      = config.getfloat('mrti','maxheat')
  fem_params['voi']          = eval(config.get('mrti','voi'))

  print 'mrti data from' , fem_params['mrti'] , 'setupfile', inisetupfile  

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


if (options.param_file != None):
  # parse the dakota input file
  fem_params = ParseInput(options.param_file)

  MatlabDriver = False
  MatlabDriver = True
  if(MatlabDriver):

    # write out for debug
    matlabInputFileName = '%s.mat' % (options.param_file)
    print matlabInputFileName 
    scipyio.savemat( matlabInputFileName , fem_params['cv'] )

    # setup any needed paths
    os.system( './analytic/dakmatlab setup workspace ' )
    matlabcommand  = './analytic/dakmatlab %s %s' %  (options.param_file,sys.argv[3])
    print matlabcommand  
    os.system( matlabcommand )
    objfunction = 0.0
  else:
    # FIXME link needed directories
    linkDirectoryList = ['occa','libocca','meshes']
    for targetDirectory in linkDirectoryList:
      linkcommand = 'ln -sf %s/%s .' % (brainNekDIR,targetDirectory )
      print linkcommand 
      os.system(linkcommand )

    # execute the rosenbrock analysis as a separate Python module
    print "Running BrainNek..."
    brainNekWrapper(**fem_params)
    
    # write objective function back to Dakota
    objfunction = ComputeObjective(**fem_params)

    print "current objective function: ",objfunction 
    fileHandle = file(sys.argv[3],'w')
    fileHandle.write('%f\n' % objfunction )
    fileHandle.flush(); fileHandle.close();

else:
  parser.print_help()
  print options
