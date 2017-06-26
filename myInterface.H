// Preprocessor
# ifndef myInterface_H
# define myInterface_H 1

// Namespace
using namespace Foam; 

// =============================================================================
//                                                            MYINTERFACE Class                                                      
// =============================================================================
//! With this abstract class it is possibile to allocate all the implemented 
//! aeroelastic interface templates as one using C++ polymorphism. A higher 
//! level class such as myInterface is also needed in order not to use only 
//! pointer allocation. Only the essental methods to be used externally are 
//! defined as virtual. 
//! TODO: More flexible handling of the interface classes with the possibility
//! of allocating multiple instances of it to be coupled with boundary patches. 
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//! 
//! \brief Abstract class to collect all the implemented aeroelastic interfaces
// =============================================================================
class MYINTERFACE
{
  
 // ----------------------------------------------------------------------------
    public:
 // ----------------------------------------------------------------------------  

    //! Destructor
    virtual ~MYINTERFACE( ) { };

 // -- Reference methods -------------------------------------------------------   
        
    //! Pointer to fluid solver
    virtual mySolver* fluid( ) = 0;
 
 // -- Advance-in-(pseudo)time operators methods -------------------------------   
        
    //! Iterate
    virtual void iterate( ) = 0;

    //! Statistics
    virtual void statistics( ) = 0;
    
    //! Operator overloading
    virtual void operator++(int) = 0;
    
 // -- Input-output methods ----------------------------------------------------     
     
    //! Rigid-Body Aerodynamic Forces 
    virtual vector& Fa( ) = 0;

    //! Rigid-Body Aerodynamic Moments
    virtual vector& Ma( ) = 0;

    //! Generalized Aerodynamic Forces
    virtual scalarField& Qa( ) = 0;    

    //! Rigid-Body Displacements
    virtual vector& s( ) = 0;

    //! Rigid-Body Orientation
    virtual vector& psi( ) = 0;   
    
    //! Rigid-Body Velocity
    virtual vector& sdot( ) = 0;

    //! Rigid-Body Orientation Velocity
    virtual vector& psidot( ) = 0;     
    
    //! Generalized Displacements
    virtual scalarField& q( ) = 0; 
    
    //! Generalized Velocities
    virtual scalarField& qdot( ) = 0;      
     
    //! Solver
    virtual word& solver( ) = 0; 
    
    //! Reference geometrical and freestream quantities
    virtual scalarField reference( ) = 0;
     
 // ----------------------------------------------------------------------------
    private:
 // ----------------------------------------------------------------------------  

};

// List of the implemented aeroelastic interfaces inheriting MYINTERFACE class
# include "myRigid.H"
# include "myModal.H"
//# include "myCodeAster.H"
//# include "myMBDyn.H"

// =============================================================================
//                                                            myInterface Class                                                      
// =============================================================================
//! This class provides the methods to read the controlDict file, select at 
//! runtime and instantiate which aeroelastic interface is to be used. This is a 
//! higher level non-virtual version of MYINTERFACE class. Only the essental 
//! methods to be used externally are wrapped. 
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//! 
//! \brief Class to select at runtime the desired aeroelastic interface
// =============================================================================
class myInterface
{

 // ----------------------------------------------------------------------------
    public:
 // ----------------------------------------------------------------------------
    
    //! Constructor
    myInterface( mySolver& fluid ) 
    {  
        // Read from file the aeroelatic interface to be set up
        _tag = "off";
        if ( fluid.time().controlDict().found("interface") ) 
        _tag = word( fluid.time().controlDict().lookup("interface") );
        
        // Dynamic memory allocation of the desired aeroelastic interface
        // Rigid
        if ( _tag == "Rigid" || _tag == "R" )
        {
            _interface = new myRigid( fluid );
        }
        // Modal
        else if ( _tag == "Modal" || _tag == "M" )
        {
            _interface = new myModal( fluid );
        }
        /*
        // MBDyn
        else if ( _tag == "MBDyn  || _tag == "MBD" )
        {
            _interface = new myMBDyn( fluid );
        }
        // CodeAster
        else if ( _tag == "CodeAster" || _tag == "CA" )
        {
            _interface = new myCodeAster( fluid );
        }*/
        else 
        {
            _tag = "off";
            _interface = NULL;
            //Info << "WARNING: No interface selected..." << endl;
        }
        
        // REMARK: Each single aeroelastic interface class will read a dedicated
        // section in a interfaceDict file in which it is possible to input all
        // the required parameters
    }
 
    //! Destructor
    ~myInterface( ) { if( _tag != "off" ) { delete _interface; } };  
    
 // -- Reference methods -------------------------------------------------------   
    
    //! Fluid
    inline mySolver* fluid( ) { if( _tag != "off" ) { return _interface->fluid(); } else { return NULL; } }  
    
 // -- Advance-in-(pseudo)time operators methods -------------------------------   
        
    //! Iterate
    inline void iterate( ) { if( _tag != "off" ) { _interface->iterate(); } }

    //! Statistics
    inline void statistics( ) { if( _tag != "off" ) { _interface->statistics(); } }
    
    //! Operator overloading
    inline void operator++(int) { if( _tag != "off" ) { _interface->iterate(); _interface->statistics(); } }

 // -- Input-output methods ----------------------------------------------------     
     
    //! Rigid-Body Aerodynamic Forces 
    inline vector& Fa( ) { return _interface->Fa(); }   

    //! Rigid-Body Aerodynamic Moments
    inline vector& Ma( ) { return _interface->Ma(); }   

    //! Generalized Aerodynamic Forces
    inline scalarField& Qa( ) { return _interface->Qa(); }   

    //! Rigid-Body Displacements
    inline vector& s( ) { return _interface->s(); }  

    //! Rigid-Body Orientation
    inline vector& psi( ) { return _interface->psi(); }     
    
    //! Rigid-Body Velocity
    inline vector& sdot( ) { return _interface->sdot(); }  

    //! Rigid-Body Orientation Velocity
    inline vector& psidot( ) { return _interface->psidot(); }  
    
    //! Generalized Displacements
    inline scalarField& q( ) { return _interface->q(); }
    
    //! Generalized Velocities
    inline scalarField& qdot( ) { return _interface->qdot(); }
    
    //! Solver
    inline word& solver( ) { return _interface->solver(); }  

    //! Reference geometrical and freestream quantities
    inline scalarField reference( ) { return _interface->reference(); }

    //! Tag
    inline word& tag( ) { return _tag; }  
      
 // ----------------------------------------------------------------------------
    private:
 // ----------------------------------------------------------------------------  
    
    //! Interface tag
    word _tag;
    
    //! Interface pointer
    MYINTERFACE* _interface;
};

# endif
