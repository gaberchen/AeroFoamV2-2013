// Preprocessor
# ifndef myImplicit_H
# define myImplicit_H 1

// Namespace
using namespace Foam; 

// =============================================================================
//                                                             myImplicit Class                                                      
// =============================================================================
//! \class myImplicit
//!
//! This class provides an interface with OpenFOAM explicit/implicit operators. 
//!
//! \author Luca Mangani, luca.mangani@gmail.com
//!         Giulio Romanelli, giulio.romanelli@gmail.com
//!
//! \brief OpenFOAM explicit/implicit operators
// =============================================================================
class MYSOLVER;

class myImplicit : public MYSOLVER
{
      
 // ----------------------------------------------------------------------------
    public:
 // ----------------------------------------------------------------------------  
           
    //! Constructor
    myImplicit( Time& inputTime, myMesh& inputMesh )
    :       
    //! Time
    _time( inputTime ),
    
    //! Mesh
    _mesh( inputMesh ),
   
    //! Thermodynamics 
    _thermodynamics( _time ),
    
    //! Navier-Stokes equations space discretization
    _NavierStokes( _time, _mesh, _thermodynamics ),
    
    //! Turbulence model equations space discretization
    _turbulence( _NavierStokes ),
    
    //! Useful variables
    _k( 0 ),
    _residual( 1.0 )  
    // TODO
    {     
        // TODO

        // Testing
        Info << "myImplicit class not implemented yet..." << endl;
        exit(-1);
    };
     
    // Destructor
    ~myImplicit( ) { };  
    
  // -- Reference access methods ------------------------------------------------
     
    //! Time
    inline Time& time() { return _time; }
 
    //! Mesh
    inline myMesh& mesh() { return _mesh; }
 
    //! Thermodynamics
    inline myThermodynamics& thermodynamics() { return _thermodynamics; }  
     
    //! NavierStokes
    inline myNavierStokes& NavierStokes() { return _NavierStokes; }
    
    //! Turbulence
    inline myTurbulence& turbulence() { return _turbulence; }
    
    //! Iteration
    inline label& iteration() { return _k; }

    //! Residual evaluation
    inline scalar residual() { return _residual; };
    
  // -- Useful methods ---------------------------------------------------------

    //! Advance solution in time
    void iterate( ); 
    
    //! Operator overloading
    void operator++(int);
    
    //! Print statistics
    void print( );
    
    //! Write statistics
    void write( );
    
    //! Print and write statistics
    void statistics( );

    // TODO
          
 // ----------------------------------------------------------------------------
    private:
 // ----------------------------------------------------------------------------
    
 // -- Reference ---------------------------------------------------------------
    
    //! Time
    Time& _time;
    
    //! Mesh
    myMesh& _mesh;
    
    //! Thermodynamics
    myThermodynamics _thermodynamics;
    
    //! NavierStokes
    myNavierStokes _NavierStokes;
           
    //! Turbulence
    myTurbulence _turbulence;       
    
 // -- Useful variables --------------------------------------------------------       

    //! Iteration number 
    label _k;

    //! Residual
    scalar _residual;

    // TODO
};

// Implementation
# include "myImplicit.C"

# endif
