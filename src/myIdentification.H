// Utilities (to be included before preprocessor definition)
# include "myIdentification.C"

// Preprocessor
# ifndef myIdentification_H
# define myIdentification_H 1

// Namespace
using namespace Foam; 

// =============================================================================
//                                                       myIdentification Class                                                      
// =============================================================================
//! Plugin class to provide all the tools necessary for the identification of
//! linear/non-linear Reduced Order Models (ROM) of the Generalized Aerodynamic
//! Forces (GAF) due to a given generalized displacement law. This can be used:
//! a) to mumerically build the state-space (non)-linear aerodynamic reduced 
//! order model and b) to numerically build the linear-ized Ham(k, M) matrix.
//! TODO: Add the option to force re-read of solution after each input is done.
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//! 
//! \brief General platform for Reduced Order Modelling (ROM) identification
// =============================================================================
class myIdentification : public MYPLUGIN
{
  
 // ----------------------------------------------------------------------------
    public:
 // ----------------------------------------------------------------------------

    //! Constructor
    myIdentification( mySolver& fluid, myInterface& structure )
    :
    //! Pointer to the fluid solver
    _fluid( &fluid ),

    //! Pointer to the structure solver
    _structure( &structure ),
    
    //! Time
    _time( _fluid->time() ),
    
    //! Mesh
    _mesh( _fluid->mesh() ),
        
    //! Thermodynamics
    _thermodynamics( _fluid->thermodynamics() ),    

    //! NavierStokes
    _NavierStokes( _fluid->NavierStokes() ),  
    
    //! Turbulence
    _turbulence( _fluid->turbulence() )
    { 
        // Read the parameters governing the identification procedure from the
        // dedicated pluginDict dictionary
        std::string parallel = ""; if ( Pstream::nProcs() > 1 ) parallel = "/..";
        IOdictionary all( IOobject( "pluginDict", _time.constant(), _time, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE ) );  
        dictionary sub;
        if ( all.found("Identification") ) 
        {
            sub = dictionary( all.subDict("Identification") ); 
        }
        else
        {
            Info << "ERROR: Identification sub-dictionary not found in pluginDict file. Aborting..." << endl;
            exit(-1);
        } 
        
        // Read from dictionary the identification strategy (linear vs. non-linear vs. Ham)
        _strategy = "linear"; if ( sub.found("strategy") ) _strategy = word( sub.lookup("strategy") ); 
       
        // Extract reference geometrical and freestream quantities from myInterface class
        scalarField reference = _structure->reference();
        _Loo = reference[0];
        _Soo = reference[1];
        _poo = reference[2];
        _Uoo = reference[3];
        _Too = reference[4];
        
        // Read input-output parameters
        _inputs.setSize(0); if ( sub.found("inputs") ) _inputs = wordList( sub.lookup("inputs") );  
        _outputs.setSize(0); if ( sub.found("outputs") ) _outputs = wordList( sub.lookup("outputs") ); 
        _i.setSize( _inputs.size(), 0.0 ); _idot.setSize( _inputs.size(), 0.0 );
        _o.setSize( _outputs.size(), 0.0 );
        
        // Read from dictionary the input time-law options
        _parameters.setSize(0); if ( sub.found("parameters") ) _parameters = scalarListList( sub.lookup("parameters") );
        if ( _parameters.size() != _inputs.size() ) 
        {
            Info << "ERROR: Prescribed time-law not matching with selected inputs. Aborting..." << endl;
            exit(-1);            
        }
        
        // Create files for input-output statistics
        if ( Pstream::myProcNo() == 0 )
        {
            // Check if the simulation is restarting or not
            std::string filenameI = _time.path() + parallel + "/Log/Identification.Inputs.log";
            std::string filenameO = _time.path() + parallel + "/Log/Identification.Outputs.log";
            FILE* fidI;
            FILE* fidO;
            if ( _time.value() == 0 ) 
            {
                fidI = fopen( &filenameI[0], "w" );
                fidO = fopen( &filenameO[0], "w" );
            }
            else
            {
                fidI = fopen( &filenameI[0], "a" );
                fidO = fopen( &filenameO[0], "a" );
            }
            fclose( fidI );
            fclose( fidO );
        }
   
        // Reset the interface solver to plugin
        _structure->solver() = "plugin";
        this->iterate();
    };

    //! Destructor
    ~myIdentification( ) { };

 // -- Reference methods -------------------------------------------------------
        
    //! Fluid
    inline mySolver* fluid( ) { return _fluid; }  

    //! Structure
    inline myInterface* structure( ) { return _structure; }  
 
 // -- Advance-in-(pseudo)time operators methods -------------------------------
        
    //! Iterate
    void iterate( );

    //! Print
    void print( );
    
    //! Write
    void write( );

    //! Statistics
    void statistics( );
    
    //! Operator overloading
    void operator++(int);
     
 // ----------------------------------------------------------------------------
    private:
 // ----------------------------------------------------------------------------  

 // -- Reference ---------------------------------------------------------------

    //! Fluid solver
    mySolver* _fluid;

    //! Structure solver
    myInterface* _structure;
    
    //! Time
    Time& _time;
    
    //! Mesh
    myMesh& _mesh;
    
    //! Thermodynamics
    myThermodynamics& _thermodynamics;
    
    //! NavierStokes
    myNavierStokes& _NavierStokes;    

    //! Turbulence
    myTurbulence& _turbulence;
 
  // -- Parameters -------------------------------------------------------------
    
    // Strategy
    word _strategy;
    
    //! Reference quantities
    scalar _poo;
    scalar _Too;
    scalar _Uoo;
    scalar _Loo;
    scalar _Soo; 
    
    // Input-output parameters
    wordList _inputs;
    wordList _outputs;
    scalarListList _parameters;  
    
    // Auxiliary storage for input-output
    scalarList _i;
    scalarList _idot;
    scalarList _o; 
};

// Implementation
# include "myIdentification.C"

# endif
