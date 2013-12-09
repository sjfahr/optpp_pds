cimport cython 

# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
#import numpy as np
#cimport numpy as np


#from libcpp.vector   cimport vector
#from libc.stdint  cimport intptr_t

# -------------------- std::string interface ----------------------- 
cdef extern from "<string>" namespace "std":
    cdef cppclass string:
        string()
        string(char *)
        char * c_str()
# -------------------- wrapper for opencl utilities ----------------------- 
cdef extern from "occa.hpp": 
    cdef void cl_list_all_devices()
# -------------------- wrapper for setupAide ----------------------- 
cdef extern from "setupAide.hpp": 
    cdef cppclass setupAide:
        setupAide(string)
        void append(string)
        string operator[](string)
# ------------------- wrapper for brain3d ----------------------- 
cdef extern from "brain3d.hpp": 
    cdef cppclass brain3d:
        brain3d(setupAide) 
        int timeStep(double)
        float dt
        void getHostTemperature(  size_t , void *)
        void setDeviceTemperature(size_t , void *)
        # http://documen.tician.de/pyopencl/misc.html#interoperability-with-other-opencl-software
        intptr_t getTemperaturePointer()
        intptr_t setTemperaturePointer(intptr_t)
## cdef extern from "mesh_base.h": 
##     cdef cppclass MeshBase:
##         Mesh (unsigned int) 
##         void read (string) 
##         void print_info () 
## # ----------------- wrapper for libMesh::GetPot ----------------------- 
## cdef extern from "getpot.h": 
##     cdef cppclass GetPot:
##         GetPot(int , char** ) 
##         void set_value "set" (char* ,char* ) # name clash w/ intrinsic ?
## # ------------- wrapper for libMesh::EquationSystems ------------------ 
## cdef extern from "parameters.h": 
##     cdef cppclass Parameters:
##         Parameters() 
##         void print_info "print" () # name clash w/ intrinsic ?
## cdef extern from "system.h": 
##     cdef cppclass System  # forward declaration
## cdef extern from "equation_systems.h": 
##     cdef cppclass EquationSystems:
##         EquationSystems(MeshBase&) 
##         void print_info () 
##         void init () 
##         System &get_system (char *) 
##         Parameters parameters 
## # ------------- wrapper for libMesh::System ------------------ 
## cdef extern from "system.h": 
##     cdef cppclass System:
##         System (EquationSystems& , string& name, unsigned int )
##         void solve ()
## # --------------- wrapper for libMesh::ExodusII_IO -------------------- 
## cdef extern from "exodusII_io.h": 
##     cdef cppclass ExodusII_IO:
##         ExodusII_IO(MeshBase &) 
##         void write_timestep        (string& , EquationSystems&, int, double)
##         void write_element_data    (EquationSystems&)
##         void write_equation_systems(string& , EquationSystems& )
## # --------------- wrapper for ITK Imaging -------------------- 
## cdef extern from "Imaging.h": 
##     cdef cppclass Imaging:
##         Imaging(GetPot &) 
##         int SetHeaderData ( int , int , int ,
## 	                    double , double , double , 
## 	                    double , double , double ) 
##         int ProjectImagingToFEMMesh (System *, PetscVec, double)
## # --------------- wrapper for Kalman -------------------- 
## cdef extern from "KalmanFilter.h": 
##     cdef cppclass KalmanFilter:
##         KalmanFilter(EquationSystems*) 
##         int Setup ( int ) 
##         int CreateIdentityMeasurementMap(    int ) 
##         int CreateMeasurementMapFromImaging( char*, int ) 
##         int StatePredict ( int ) 
##         int CovariancePredict ( ) 
##         int StateUpdate ( char* ) 
##         int CovarianceUpdate ( ) 
##         int MeasurementMatrixTranspose( char*,double )
##         int ProjectMeasurementMatrix(PetscVec,PetscVec *)
##         int ExtractMeasurementData ( PetscVec ,char*)
##         int ExtractVarianceForPlotting ( char* ) 
##         int ExtractCoVarianceForPlotting ( char*,int ) 
##         int CreateROIAverageMeasurementMap(int,int,int,
##                                            int,int,int,int,double,double,double,
##                                                            double,double,double)
##         int CreateROINodeSetMeasurementMap(int)
##         void printSelf () 
##     cdef cppclass DenseKalmanFilter:
##         DenseKalmanFilter(EquationSystems*) 
##     cdef cppclass SparseKalmanFilter:
##         SparseKalmanFilter(EquationSystems*) 
##     cdef cppclass UncorrelatedKalmanFilter:
##         UncorrelatedKalmanFilter(EquationSystems*) 
## # --------------- c-style wrappings for interface -------------------- 
## cdef extern from "tttkUtilities.h": 
##     cdef int GenerateStructuredGrid(Mesh &,
##                                     int , int , int ,
## 				    double , double , 
## 				    double , double , 
## 				    double , double , 
##                                     vector[int] &)
##     cdef int SetupUnStructuredGrid(Mesh *,char *, int,
##                                    double, double, double, 
##                                    double, double, double, 
##                                    double, double, double, 
##                                    double, double, double)
##     cdef int SetParameterIni( Parameters *, GetPot *)
##     cdef int PetscPopErrorHandlerAndDebug( )
##     cdef int FlushStdCoutCerr( )
##     cdef System*  GetSystem (            char *, EquationSystems *)
##     cdef System*  AddExplicitSystem (    char *, EquationSystems *)
##     cdef System*  AddBackgroundSystem   (char *, EquationSystems *)
##     cdef System*  AddPennesDeltaPSystem (char *, EquationSystems *,double)
##     cdef System*  AddPennesSDASystem    (char *, EquationSystems *,double)
##     cdef System*  AddPennesRFSystem     (char *, EquationSystems *,double)
##     cdef unsigned int AddConstantMonomialVariable( System*, char*)
##     cdef unsigned int AddFirstLagrangeVariable(    System*, char*)
##     cdef int StoreSystemTimeStep( System*, int)
##     cdef int  GetSolutionVector(   System&, PetscVec *)
##     cdef int  SetSolutionVector(   System&, PetscVec  )
##     cdef int CopySolutionVector(   System&, System&   )
##     cdef int AddStorageVectors(   System*, char * , int )
##     cdef double WEIGHTEDL2Norm(   System*, char * ,System*, char * ,System*,char *)
##     cdef PetscMat GETFactorMat( PetscMat, char *,char *)
##     cdef int PetscFEMSystemUpdateTimeStep(        System *, int)
##     cdef int PetscFEMSystemSetupInitialConditions(System *)
##     cdef int PetscFEMSystemGetSolnSubVector(      System *, int,PetscVec *)
##     cdef int PennesSDASystemUpdateLaserPosition(System* ,
##                       double,double,double,double,double,double ) 
##     cdef int PennesSDASystemUpdateLaserPower(System* , double,int) 
##     cdef int PetscFEMSystemCreateNodeSetFromMask( System* ,double,int)
