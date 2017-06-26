// Preprocessor
# ifndef myTurbulence_H
# define myTurbulence_H 1

// Global variables
# define MG_TURON 1 // Turbulence model equations are embedded in Multi-Grid (0, 1)    

// Namespace
using namespace Foam; 

// =============================================================================
//                                                           MYTURBULENCE Class                                                      
// =============================================================================
//! With this abstract class it is possibile to allocate all the implemented 
//! turbulence models as one using C++ polymorphism. A higher level class such 
//! as myTurbulence is also needed in order not to use only pointer allocation.
//! Only the essental methods to be used externally (e.g. in mySolver) are 
//! defined as virtual. 
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//! 
//! \brief Abstract class to collect all the implemented turbulence models
// =============================================================================
class MYTURBULENCE
{
  
 // ----------------------------------------------------------------------------
    public:
 // ----------------------------------------------------------------------------  
    
    //! Destructor
    virtual ~MYTURBULENCE( ) { };

 // -- Reference methods -------------------------------------------------------   
        
    //! Modelling tag
    virtual inline word& tag() = 0;

    //! Wall distance
    virtual wallDist& wallDistance() = 0;
 
 // -- Space discretization operators methods ----------------------------------    
        
    //! Inviscid advection 
    virtual void advection( ) = 0;
    
    //! Laminar and turbulent diffusion
    virtual void diffusion( ) = 0;
     
    //! Source terms
    virtual void source( bool unsteady ) = 0; 
    
    //! External source terms (body forces)
    virtual void body( bool unsteady ) = 0; 
    
    //! Set to zero rhs
    virtual void resetRhs( ) = 0;
    
    //! Set to zero body rhs
    virtual void resetBody( ) = 0;

    //! Smooth rhs
    virtual void smoothRhs( label iterations, scalar epsilon ) = 0;
    
 // -- Time discretization operators methods -----------------------------------     
    
    //! Solve turbulence model  
    virtual void solve( scalar alpha, label iterations, scalar epsilon ) = 0;

    //! Store the solution at timestep (k) as (k - 1)
    virtual void store() = 0;

    //! Update the auxiliary variables and the interface between myNavierStokes and 
    //! myTurbulence classes with an update of muTur and kTur variables
    virtual void update() = 0;

    //! Activate wall functions
    virtual void wallFunctions() = 0;
    
 // -- Residuals access methods ------------------------------------------------     
    
    //! Residual
    virtual scalar residual() = 0;
  
    //! Reset residual
    virtual void resetResidual() = 0; 
     
    //! Update residual
    virtual void updateResidual( word normalization ) = 0;
    
 // -- Dual TimeStepping (DTS) methods -----------------------------------------   
  
    //! Store 1-st half and update 2-nd half of source term for DTS
    virtual void buildDTS( label half ) = 0;

 // -- Opaque containers -------------------------------------------------------
    
    //! Number of equations
    virtual label& size() = 0;
    
    //! Conservative variables at timestep k
    virtual scalarField& conservative( label ic ) = 0;

    //! Conservative variables at timestep k-1
    virtual scalarField& conservative_o( label ic ) = 0;   

    //! Body forces
    virtual scalarField& body( label ic ) = 0;   
    
    //! Rhs array
    virtual scalarField& rhs( label ic ) = 0;   
    
    //! Correct boundary conditions
    virtual void correctBoundaryConditions( ) = 0;        
     
 // ----------------------------------------------------------------------------
    private:
 // ----------------------------------------------------------------------------  

};

// List of the implemented turbulence models inheriting MYTURBULENCE class
# include "mySpalartAllmaras.H"
# include "myKappaOmega.H"

// =============================================================================
//                                                           myTurbulence Class                                                      
// =============================================================================
//! This class provides the methods to read the controlDict file, select at 
//! runtime and instantiate which turbulence model is to be used. This is a 
//! higher level non-virtual version of MYTURBULENCE class. Only the essental
//! methods to be used externally (e.g. in mySolver) are wrapped. In this way 
//! nothing is allocated if tag is set to off.
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//! 
//! \brief Class to select at runtime the desired turbulence model
// =============================================================================
class myTurbulence
{

 // ----------------------------------------------------------------------------
    public:
 // ----------------------------------------------------------------------------
    
    //! Constructor
    myTurbulence( myNavierStokes& inputNavierStokes ) 
    {  
        // Read from file the turbulence model (vs off)
        _tag = "off";
        if ( inputNavierStokes.time().controlDict().found("turbulence") ) 
        _tag = word( inputNavierStokes.time().controlDict().lookup("turbulence") );
        
        // Turbulence model is not active on coarse mesh levels 
        // (it should be disabled by default, e.g. MG_TURON = 1)
        # if MG_TURON == 0
        if ( inputNavierStokes.mesh().tag() != "*" ) _tag = "off";
        # endif 
        
        // Dynamic memory allocation of the desired turbulence model
        // Spalart-Allmaras 1-eq turbulence model
        if ( _tag == "SpalartAllmaras" || _tag == "SA" )
        {
            _turbulence = new mySpalartAllmaras( inputNavierStokes );
            _size = 1;
        }
        // Kappa-Omega 2-eqs turbulence model 
        else if ( _tag == "KappaOmega" || _tag == "KW" )
        {
            _turbulence = new myKappaOmega( inputNavierStokes );
            _size = 2;
        }
        else 
        {
            _tag = "off";
            _size = 0;
            _dummy = scalarField( 1, 0.0 );
            _turbulence = NULL;
        }
    }
 
    //! Destructor
    ~myTurbulence( ) { if ( _tag != "off" ) delete _turbulence; };  
    
 // -- Reference methods -------------------------------------------------------   
        
    //! Modelling tag (redundant also in each turbulence model implementation)
    inline word& tag() { return _tag; }  

 // -- Space discretization operators methods ----------------------------------    

    //! Update the wall distance field y for ALE formulation with changing meshes 
    inline void updateWallDistance() { if ( _tag != "off" ) _turbulence->wallDistance().correct(); }
        
    //! Inviscid advection 
    inline void advection() { if ( _tag != "off" ) _turbulence->advection(); }
    
    //! Laminar and turbulent diffusion
    inline void diffusion() { if ( _tag != "off" ) _turbulence->diffusion(); }
     
    //! Source terms
    inline void source( bool unsteady = false ) { if ( _tag != "off" ) _turbulence->source( unsteady ); }
    
    //! External source terms (body forces)
    inline void body( bool unsteady = false ) { if ( _tag != "off" ) _turbulence->body( unsteady ); } 
    
    //! Set to zero rhs
    inline void resetRhs() { if ( _tag != "off" ) _turbulence->resetRhs( ); } 
   
    //! Set to zero body rhs
    inline void resetBody() { if ( _tag != "off" ) _turbulence->resetBody( ); } 

    //! Smooth rhs
    inline void smoothRhs( label iterations, scalar epsilon ) { if ( _tag != "off" ) _turbulence->smoothRhs( iterations, epsilon ); } 
    
 // -- Time discretization operators methods -----------------------------------      
    
    //! Solve turbulence model  
    inline void solve( scalar alpha, label iterations, scalar epsilon ) { if ( _tag != "off" ) _turbulence->solve( alpha, iterations, epsilon ); }
    
    //! Store the solution at timestep (k) as (k - 1)
    inline void store() { if ( _tag != "off" ) _turbulence->store(); }    
    
    //! Update the auxiliary variables and the interface between myNavierStokes and 
    //! myTurbulence classes with an update of muTur and kTur variables
    inline void update() { if ( _tag != "off" ) _turbulence->update(); } 

    //! Activate wall functions
    inline void wallFunctions(){ if ( _tag != "off" ) _turbulence->wallFunctions(); } 
    
 // -- Residuals access methods ------------------------------------------------     
    
    //! Residual (the check on the tag must be done externally)
    inline scalar residual() { if ( _tag != "off" ) return _turbulence->residual(); else return -1.0; } 
  
    //! Reset residual
    inline void resetResidual() { if ( _tag != "off" ) _turbulence->resetResidual(); } 
     
    //! Update residual
    inline void updateResidual( word normalization ) { if ( _tag != "off" ) _turbulence->updateResidual( normalization ); } 
    
 // -- Dual TimeStepping (DTS) methods -----------------------------------------   
  
    //! Store 1-st half and update 2-nd half of source term for DTS
    inline void buildDTS( label half ) { if ( _tag != "off" ) _turbulence->buildDTS( half ); } 

 // -- Opaque containers -------------------------------------------------------
    
    //! Number of equations
    inline label& size() { return _size; }

    //! Conservative variables at timestep k 
    inline scalarField& conservative( label ic ) { if ( _tag != "off" ) return _turbulence->conservative( ic ); else return _dummy; } 

    //! Conservative variables at timestep k-1
    inline scalarField& conservative_o( label ic ) { if ( _tag != "off" ) return _turbulence->conservative_o( ic ); else return _dummy; }    

    //! Body forces
    inline scalarField& body( label ic ) { if ( _tag != "off" ) return _turbulence->body( ic ); else return _dummy; } 

    //! Rhs array
    inline scalarField& rhs( label ic ) { if ( _tag != "off" ) return _turbulence->rhs( ic ); else return _dummy; } 
   
    //! Correct boundary conditions
    void correctBoundaryConditions( ) { if ( _tag != "off" ) _turbulence->correctBoundaryConditions( ); }
     
 // ----------------------------------------------------------------------------
    private:
 // ----------------------------------------------------------------------------  

    //! Modelling tag (redundant also in each turbulence model implementation)
    word _tag;
    label _size;
    scalarField _dummy;
    
    //! Turbulence model pointer
    MYTURBULENCE* _turbulence;
};

# endif
