// Preprocessor
# ifndef mySolver_H
# define mySolver_H 1

// Namespace
using namespace Foam; 

// =============================================================================
//                                                               MYSOLVER Class                                                      
// =============================================================================
//! With this abstract class it is possibile to allocate all the implemented 
//! solver templates as one using C++ polymorphism. A higher level class such 
//! as mySolver is also needed in order not to use only pointer allocation.
//! Only the essental methods to be used externally are defined as virtual. 
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//! 
//! \brief Abstract class to collect all the implemented solvers
// =============================================================================
class MYSOLVER
{
  
 // ----------------------------------------------------------------------------
    public:
 // ----------------------------------------------------------------------------  
    
    //! Destructor
    virtual ~MYSOLVER( ) { };

 // -- Reference methods -------------------------------------------------------   
        
    //!  Time
    virtual Time& time( ) = 0;
    
    //! Mesh
    virtual myMesh& mesh() = 0;
    
    //! Thermodynamics
    virtual myThermodynamics& thermodynamics() = 0;
     
    //! NavierStokes
    virtual myNavierStokes& NavierStokes() = 0;
    
    //! Turbulence
    virtual myTurbulence& turbulence() = 0;
    
    //! Iteration
    virtual label& iteration() = 0; 

    //! Iteration
    virtual scalar residual() = 0; 
     
 // -- Advance-in-(pseudo)time operators methods -------------------------------   
        
    //! Iterate
    virtual void iterate( ) = 0;

    //! Statistics
    virtual void statistics( ) = 0;
    
    //! Operator overloading
    virtual void operator++(int) = 0;
     
 // ----------------------------------------------------------------------------
    private:
 // ----------------------------------------------------------------------------  

};

// List of the implemented solver templates inheriting MYSOLVER class
# include "myTimeStepping.H"
//# include "myMultiGrid.H"
# include "myImplicit.H"

// =============================================================================
//                                                               mySolver Class                                                      
// =============================================================================
//! This class provides the methods to read the controlDict file, select at 
//! runtime and instantiate which solver template is to be used. This is a 
//! higher level non-virtual version of MYSOLVER class. Only the essental 
//! methods to be used externally are wrapped. 
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//! 
//! \brief Class to select at runtime the desired solver
// =============================================================================
class mySolver
{

 // ----------------------------------------------------------------------------
    public:
 // ----------------------------------------------------------------------------
    
    //! Constructor
    mySolver( Time& inputTime, myMesh& inputMesh ) 
    {  
        // Read from file the solver
        _tag = "default";
        if ( inputTime.controlDict().found("solver") ) 
        _tag = word( inputTime.controlDict().lookup("solver") );
        
        // Read from file if Log folder must be cleared
        word clear = "off";
        if ( inputTime.controlDict().found("clear") ) 
        clear = word( inputTime.controlDict().lookup("clear") );
        std::string parallel = ""; if ( Pstream::nProcs() > 1 ) parallel = "/..";
        std::string filename = inputTime.path() + parallel + "/Log/*";
        if ( clear == "all" ) system( "rm " + filename + " > /dev/null 2>&1" );
                
        // Dynamic memory allocation of the desired solver
        // Multi-Grid
        //if ( _tag == "MultiGrid" || _tag == "MG" )
        //{
        //    _solver = new myMultiGrid( inputTime, inputMesh );
        //}
        // Implicit
        //else if ( _tag == "Implicit" || _tag == "I" )
        //{
        //    _solver = new myImplicit( inputTime, inputMesh );
        //}
        // Time-Stepping
        //else 
        //{
            _solver = new myTimeStepping( inputTime, inputMesh );
        //}

    }
 
    //! Destructor
    ~mySolver( ) { delete _solver; };  
    
 // -- Reference methods -------------------------------------------------------   
    
    //!  Time
    inline Time& time() { return _solver->time(); }
    
    //! Mesh
    inline myMesh& mesh() { return _solver->mesh(); }
    
    //! Thermodynamics
    inline myThermodynamics& thermodynamics() { return _solver->thermodynamics(); }
     
    //! NavierStokes
    inline myNavierStokes& NavierStokes() { return _solver->NavierStokes(); }
    
    //! Turbulence
    inline myTurbulence& turbulence() { return _solver->turbulence(); }   

    //! Iteration
    inline label& iteration() { return _solver->iteration(); }
    
    //! Residual
    inline scalar residual() { return _solver->residual(); }
        
 // -- Advance-in-(pseudo)time operators methods -------------------------------   
        
    //! Iterate
    inline void iterate() { _solver->iterate(); }

    //! Statistics
    inline void statistics() { _solver->statistics(); }
    
    //! Operator overloading
    inline void operator++(int) { _solver->iterate(); _solver->statistics(); }
    
 // ----------------------------------------------------------------------------
    private:
 // ----------------------------------------------------------------------------  
    
    //! Solver tag
    word _tag;
    
    //! Solver pointer
    MYSOLVER* _solver;
};

# endif
