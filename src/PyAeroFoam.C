// System include
# include "Python.h"

// Project include
# include "myOpenFOAM.H"
# include "myMesh.H"
# include "myThermodynamics.H"
# include "myNavierStokes.H"
# include "myTurbulence.H"
# include "mySolver.H"
# include "myInterface.H"
# include "myPlugin.H"

// =============================================================================
//                                                     PyAeroFoam Wrapper (1/2)                                                 
// =============================================================================
//! \class PyAeroFoam
//!
//! Utilities to build a Python wrapper for AeroFoam. More in particular it is 
//! possibile to find: a) global variables definitions, b) pre-wrapper functions
//! to simiplify the use of global variables and c) Python wrapper functions.
//! REMARK: The prefix is always written as "Py" to prevent confusion between
//!         e.g. myClass and its wrapper PyClass.
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//!
//! \brief Utilities to build a Python wrapper for AeroFoam 
// =============================================================================

//! Global variables definitions as pointers to the following global types: 
//! (built-in time and mesh), myMesh, mySolver, myInterface, myPlugin...
Time* PyTime;
fvMesh* PyRead;
myMesh* PyMesh;
mySolver* PySolver;
myInterface* PyInterface;
myPlugin* PyPlugin;

//! Definition of the id of each global class
int IDTIME;
int IDMESH;
int IDSOLV;
int IDINTE;
int IDPLUG;

// =============================================================================
//                                                                   PyTime_New                                                
// =============================================================================
//! Pre-wrapper function to allocate OpenFOAM built-in class Time
Time* PyTime_New( char* rootPath, char* caseName )
{
    Time* memory = new Time( word("controlDict"), fileName(rootPath), fileName(caseName) );
    return memory;
}

// =============================================================================
//                                                                PyTime_Delete                                              
// =============================================================================
//! Pre-wrapper function to deallocate OpenFOAM built-in class Time
void PyTime_Delete( Time* memory )
{
    delete memory;
}

// =============================================================================
//                                                                   PyRead_New                                                
// =============================================================================
//! Pre-wrapper function to allocate OpenFOAM built-in class fvMesh
fvMesh* PyRead_New( Time& time )
{
    fvMesh* memory = new fvMesh( IOobject( fvMesh::defaultRegion, time.timeName(), time, IOobject::MUST_READ ) );  
    return memory;
}

// =============================================================================
//                                                                PyRead_Delete                                              
// =============================================================================
//! Pre-wrapper function to deallocate OpenFOAM built-in class fvMesh
void PyRead_Delete( fvMesh* memory )
{
    delete memory;
}

// =============================================================================
//                                                                   PyMesh_New                                                
// =============================================================================
//! Pre-wrapper function to allocate class myMesh
myMesh* PyMesh_New( Time& time, fvMesh& read )
{
    myMesh* memory = new myMesh( time, read );  
    return memory;
}

// =============================================================================
//                                                                PyMesh_Delete                                              
// =============================================================================
//! Pre-wrapper function to deallocate class myMesh
void PyMesh_Delete( myMesh* memory )
{
    delete memory;
}

// =============================================================================
//                                                                 PySolver_New                                                
// =============================================================================
//! Pre-wrapper function to allocate class mySolver
mySolver* PySolver_New( Time& time, myMesh& mesh )
{
    mySolver* memory = new mySolver( time, mesh );  
    return memory;
}

// =============================================================================
//                                                              PySolver_Delete                                              
// =============================================================================
//! Pre-wrapper function to deallocate class mySolver
void PySolver_Delete( mySolver* memory )
{
    delete memory;
}

// =============================================================================
//                                                              PyInterface_New                                                
// =============================================================================
//! Pre-wrapper function to allocate class myInterface
myInterface* PyInterface_New( mySolver& solver )
{
    myInterface* memory = new myInterface( solver );  
    return memory;
}

// =============================================================================
//                                                           PyInterface_Delete                                              
// =============================================================================
//! Pre-wrapper function to deallocate class myInterface
void PyInterface_Delete( myInterface* memory )
{
    delete memory;
}

// =============================================================================
//                                                                 PyPlugin_New                                                
// =============================================================================
//! Pre-wrapper function to allocate class myPlugin
myPlugin* PyPlugin_New( mySolver& solver, myInterface& interface )
{
    myPlugin* memory = new myPlugin( solver, interface );  
    return memory;
}

// =============================================================================
//                                                              PyPlugin_Delete                                              
// =============================================================================
//! Pre-wrapper function to deallocate class myPlugin
void PyPlugin_Delete( myPlugin* memory )
{
    delete memory;
}

// =============================================================================
//                                                           PyWrapper_OpenFOAM                                           
// =============================================================================
//! Set the root path and case name. In this way it is possibile to run a case 
//! from differents folders, e.g. to set-up a parametric case queue.
PyObject* PyWrapper_OpenFOAM( PyObject *self, PyObject *args ) 
{
    // Variables definition
    char *rootPath, *caseName;
    
    // Parser and set global id
    int id = -1;
    PyArg_ParseTuple( args, "iss;", &id, &rootPath, &caseName );
    IDTIME = id;
    
    // Memory allocation 
    PyTime = PyTime_New( rootPath, caseName );
    PyRead = PyRead_New( *PyTime );
    
    // Return
    return Py_None;
} 

// =============================================================================
//                                                               PyWrapper_Mesh                                              
// =============================================================================
//! Wrapper function for myMesh dynamic allocation
PyObject* PyWrapper_Mesh( PyObject *self, PyObject *args ) 
{
    // Parser and set global id
    int id = -1;
    PyArg_ParseTuple( args, "i;", &id );
    IDMESH = id;

    // Memory allocation 
    PyMesh = PyMesh_New( *PyTime, *PyRead );
    
    // Return
    return Py_None;
} 

// =============================================================================
//                                                             PyWrapper_Solver                                          
// =============================================================================
//! Wrapper function for mySolver dynamic allocation
PyObject* PyWrapper_Solver( PyObject *self, PyObject *args ) 
{
    // Parser and set global id
    int id = -1;
    PyArg_ParseTuple( args, "i;", &id );
    IDSOLV = id;

    // Memory allocation 
    PySolver = PySolver_New( *PyTime, *PyMesh );

    // Return
    return Py_None;  
} 

// =============================================================================
//                                                          PyWrapper_Interface                                            
// =============================================================================
//! Wrapper function for myInterface dynamic allocation
PyObject* PyWrapper_Interface( PyObject *self, PyObject *args ) 
{
    // Parser and set global id
    int id = -1;
    PyArg_ParseTuple( args, "i;", &id );
    IDINTE = id;

    // Memory allocation 
    PyInterface = PyInterface_New( *PySolver );

    // Return
    return Py_None;  
} 

// =============================================================================
//                                                            PyWrapper_Plugin                                              
// =============================================================================
//! Wrapper function for myPlugin dynamic allocation
PyObject* PyWrapper_Plugin( PyObject *self, PyObject *args ) 
{
    // Parser and set global id
    int id = -1;
    PyArg_ParseTuple( args, "i;", &id );
    IDPLUG = id;

    // Memory allocation 
    PyPlugin = PyPlugin_New( *PySolver, *PyInterface );

    // Return
    return Py_None;  
}

// =============================================================================
//                                                               PyWrapper_Next                                            
// =============================================================================
//! Wrapper function to advance in time the id-th class
PyObject* PyWrapper_Next( PyObject *self, PyObject *args ) 
{
    // Parser and set global id
    int id = -1;
    PyArg_ParseTuple( args, "i;", &id );
    
    //--------------------------------------------------------------------------
    // Time
    //--------------------------------------------------------------------------
    if ( id == IDTIME ) 
    {
        PyTime->operator++(1);
        PyTime->write();    
        if ( PyTime->end() ) 
        { 
            return Py_False;
        }
        else
        {
            return Py_True;
        }    
    }
    //--------------------------------------------------------------------------
    // Mesh
    //--------------------------------------------------------------------------
    else if ( id == IDMESH )
    {
        // Do nothing... 
        PyMesh->operator++(1);
        return Py_None;    
    }
    //--------------------------------------------------------------------------
    // Solver
    //--------------------------------------------------------------------------
    else if ( id == IDSOLV )
    {
        PySolver->operator++(1);
        return Py_None;
    }
    //--------------------------------------------------------------------------
    // Interface
    //--------------------------------------------------------------------------
    else if ( id == IDINTE )
    {
        PyInterface->operator++(1);
        return Py_None;
    }  
    //--------------------------------------------------------------------------
    // Plugin
    //--------------------------------------------------------------------------
    else if ( id == IDPLUG )
    {
        PyPlugin->operator++(1);
        return Py_None; 
    } 
    //--------------------------------------------------------------------------
    // Exceptions
    //--------------------------------------------------------------------------
    else
    {
        // Do nothing... 
        return Py_None;  
    }      
}

// =============================================================================
//                                                               PyWrapper_Free                                            
// =============================================================================
//! Wrapper function to free allocated memory
PyObject* PyWrapper_Free( PyObject *self, PyObject *args ) 
{
    // Parser and set global id
    int id = -1;
    PyArg_ParseTuple( args, "i;", &id );
    
    //--------------------------------------------------------------------------
    // Time
    //--------------------------------------------------------------------------
    if ( id == IDTIME ) 
    {   
        delete PyRead;
        delete PyTime;
        return Py_None;  
    }
    //--------------------------------------------------------------------------
    // Mesh
    //--------------------------------------------------------------------------
    else if ( id == IDMESH )
    {
        delete PyMesh;
        return Py_None;    
    }
    //--------------------------------------------------------------------------
    // Solver
    //--------------------------------------------------------------------------
    else if ( id == IDSOLV )
    {
        delete PySolver;
        return Py_None;    
    }
    //--------------------------------------------------------------------------
    // Interface
    //--------------------------------------------------------------------------
    else if ( id == IDINTE )
    {
        delete PyInterface;
        return Py_None; 
    }  
    //--------------------------------------------------------------------------
    // Plugin
    //--------------------------------------------------------------------------
    else if ( id == IDPLUG )
    {
        delete PyPlugin;
        return Py_None; 
    } 
    //--------------------------------------------------------------------------
    // Exceptions
    //--------------------------------------------------------------------------
    else
    {
        // Do nothing... 
        return Py_None;  
    }      
}

// =============================================================================
//                                                            PyWrapper_Methods                                              
// =============================================================================
//! Python methods to be wrapped in the module
static PyMethodDef PyWrapper_Methods[] = { { "OpenFOAM",  PyWrapper_OpenFOAM,  METH_VARARGS, NULL },
                                           { "Mesh",      PyWrapper_Mesh,      METH_VARARGS, NULL },  
                                           { "Solver",    PyWrapper_Solver,    METH_VARARGS, NULL }, 
                                           { "Interface", PyWrapper_Interface, METH_VARARGS, NULL }, 
                                           { "Plugin",    PyWrapper_Plugin,    METH_VARARGS, NULL }, 
                                           { "Next",      PyWrapper_Next,      METH_VARARGS, NULL }, 
                                           { "Free",      PyWrapper_Free,      METH_VARARGS, NULL }, 
                                           { NULL, NULL } };

// =============================================================================
//                                                               initPyAeroFoam                                              
// =============================================================================
//! Initialization function called when import <module> command is issued
extern "C" void initPyAeroFoam() 
{
    // Initialization
    PyObject *self;
    self = Py_InitModule("PyAeroFoam", PyWrapper_Methods);
}    
