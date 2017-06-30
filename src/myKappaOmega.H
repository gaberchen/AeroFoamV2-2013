// Preprocessor
# ifndef myKappaOmega_H
# define myKappaOmega_H 1

// High resolution advection fluxes
# define KW_HIRE 0      // Activate high resolution fluxes (0, 1)
# define KW_SMALL 1e-10 // Minium value for turbulent kappa and omega
# define KW_LINFIX 0.10 // Entropy fix on linear field (0, 0.1)

// Namespace
using namespace Foam; 

// =============================================================================
//                                                     myKappaOmegaPatch Struct
// =============================================================================
//! \class myKappaOmegaPatch
//!
//! This struct provides a blackboard to copy data on for cyclic and processor
//! boundary patches and therefore to make communication easier and seamless. 
//! To be used in conjunction with the class myKappaOmega.
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//!
//! \brief Blackboard to copy data on for cyclic and processor boundary patches
// =============================================================================
typedef struct{ scalarField rho_R;  
                vectorField U_R;  
                scalarField kappa_R; scalarField kappa_RR;
                scalarField omega_R; scalarField omega_RR;
              } myKappaOmegaPatch;

// =============================================================================
//                                                           myKappaOmega Class                                                      
// =============================================================================
//! \class myKappaOmega
//!
//! This class provides the space discretization operators for the turbulence 
//! model equations, categorized into advection, diffusion and source access
//! methods. To enable a high resolution discretization set the flag KW_HIRE 
//! to 1 by means of suitable pre-processor call.
//! This class is specialized for the Kappa-Omega SST 1-eq turbulence model 
//! [Fluent User's Guide], but its structure is general to any RANS turbulence 
//! model.
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//!
//! \brief Space discretization operators for the turbulence model equations
// =============================================================================
class MYTURBULENCE;

class myKappaOmega : public MYTURBULENCE
{
  
 // ----------------------------------------------------------------------------
    public:
 // ----------------------------------------------------------------------------  
 
    //! Constructor
    myKappaOmega( myNavierStokes& inputNavierStokes )
    :
    //! Tag
    _tag( "KappaOmega" ),
    _size( 2 ),
     
    //! Reference
    _NavierStokes( inputNavierStokes ),
    _time( inputNavierStokes.time() ),
    _mesh( inputNavierStokes.mesh() ),
    _thermodynamics( inputNavierStokes.thermodynamics() ),
    
    //! kappa and omega arrays
    _kappa( IOobject("kappa", _time.timeName(), _mesh.mesh(), IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE ), _mesh.mesh() ),
    _kappa_o( IOobject("kappa_o", _time.timeName(), _mesh.mesh(), IOobject::READ_IF_PRESENT, IOobject::NO_WRITE ), _kappa, _kappa.boundaryField().types() ),
    _gradKappa( fvc::grad( _kappa ) ),
    _omega( IOobject("omega", _time.timeName(), _mesh.mesh(), IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE ), _mesh.mesh() ),
    _omega_o( IOobject("omega_o", _time.timeName(), _mesh.mesh(), IOobject::READ_IF_PRESENT, IOobject::NO_WRITE ), _omega, _omega.boundaryField().types() ),
    _gradOmega( fvc::grad( _omega ) ),
    
    //! Residuals
    _residualKappa( -1.0 ),
    _residualOmega( -1.0 ),
    _maxResidualKappa( 1.0e-16 ),
    _maxResidualOmega( 1.0e-16 ),
    
    //! Rhs and lhs arrays
    _rhsKappa( _mesh.V().size(), 0.0 ),
    _rhsOmega( _mesh.V().size(), 0.0 ),
    _lhsKappa( _mesh.V().size(), 0.0 ),
    _lhsOmega( _mesh.V().size(), 0.0 ),
        
    // Wall distances
    _d( _mesh.mesh() ),
    
    //! Constants
    _k( 0.41 ),
    _alphaKappa1( 0.85034 ),
    _alphaKappa2( 1.0 ),    
    _alphaOmega1( 0.5 ),
    _alphaOmega2( 0.85616 ), 
    _gamma1( 0.5532 ),
    _gamma2( 0.4403 ),
    _beta1( 0.075 ),
    _beta2( 0.0828 ),
    _betaStar( 0.09 ),    
    _a1( 0.31 ),
    _c1( 10.0 ),
    _C( 5.5 ),
    _E( Foam::exp( _C*_k ) ),
            
    //! Dual TimeStepping
    _dtsKappa( _mesh.V().size(), 0.0 ),
    _dtsOmega( _mesh.V().size(), 0.0 ),
    
    //! External source terms (body forces)
    _bodyKappa( _mesh.V().size(), 0.0 ),
    _bodyOmega( _mesh.V().size(), 0.0 )
    {  
        // Read from file, if present, the turbulence model (SpalartAllmaras vs off)
        if ( _time.controlDict().found("turbulence") ) 
        _tag = word( _time.controlDict().lookup("turbulence") );
        
        // Check bounds
        _kappa.internalField() = max( _kappa.internalField(), KW_SMALL );
        _omega.internalField() = max( _omega.internalField(), KW_SMALL );
        
        
    };
     
    //! Destructor
    ~myKappaOmega( ) { };  

 // -- Reference access methods ------------------------------------------------
 
    //! Modelling tag
    inline word& tag() { return _tag; }

    //! NavierStokes
    inline myNavierStokes& NavierStokes() { return _NavierStokes; }
    
    //! Time
    inline Time& time() { return _time; }
 
    //! Mesh
    inline myMesh& mesh() { return _mesh; }
 
    //! Thermodynamics
    inline myThermodynamics& thermodynamics() { return _thermodynamics; }

    //! Wall distance 
    inline wallDist& wallDistance() { return _d; }

 // -- Space discretization operators methods ----------------------------------

    //! Inviscid advection 
    void advection( );
    
    //! Laminar and turbulent diffusion
    void diffusion( );
     
    //! Source terms
    void source( bool unsteady ); 
    
    //! External source terms (body forces)
    void body( bool unsteady );
    
 // -- Kappa-Omega conservative variables access methods -----------------------
   
    //! kappa at timestep (k)
    inline volScalarField& kappa() { return _kappa; }
    
    //! kappa at timestep (k - 1)
    inline volScalarField& kappa_o() { return _kappa_o; }
    
    //! omega at timestep (k)
    inline volScalarField& omega() { return _omega; }
    
    //! omega at timestep (k - 1)
    inline volScalarField& omega_o() { return _omega_o; } 
       
    //! Turbulent viscosity
    volScalarField muTur();
    
    //! Turbulent kinetic energy
    volScalarField kTur();
     
    //! Solve K-W turbulence model  
    void solve( scalar alpha, label iterations, scalar epsilon );

    //! Store the solution at timestep (k) as (k - 1)
    void store();    

    //! Update the auxiliary variables and the interface between myNavierStokes and 
    //! myKappaOmega classes with an update of muTur and kTur variables
    void update();

    //! Activate wall functions
    void wallFunctions();    
   
 // -- Residuals access methods ------------------------------------------------     
    
    //! Residual
    inline scalar residual() { return max( _residualKappa, _residualOmega ); }
  
    //! Reset residual
    void resetResidual(); 
     
    //! Update residual
    void updateResidual( word normalization );
  
 // -- Rhs and lhs access methods ---------------------------------------------- 
  
    //! Rhs
    inline scalarField& rhsKappa() { return _rhsKappa; }
    inline scalarField& rhsOmega() { return _rhsOmega; }
    
    //! Lhs
    inline scalarField& lhsKappa() { return _lhsKappa; }
    inline scalarField& lhsOmega() { return _lhsOmega; }
    
    //! Set to zero rhs
    void resetRhs();
    
    //! Smooth rhs
    void smoothRhs( label iterations, scalar epsilon );
    
    //! Body rhs
    inline scalarField& bodyKappa() { return _bodyKappa; }
    inline scalarField& bodyOmega() { return _bodyOmega; }
    
    //! Set to zero body rhs
    void resetBody();
    
 // -- Kappa-Omega constants ---------------------------------------------------
  
    //! Von Karman k
    inline scalar& k() { return _k; }
    
    //! alphaKappa1
    inline scalar alphaKappa1() { return _alphaKappa1; }
    
    //! alphaKappa2
    inline scalar alphaKappa2() { return _alphaKappa2; }
    
    //! alphaOmega1
    inline scalar alphaOmega1() { return _alphaOmega1; }
    
    //! alphaOmega2
    inline scalar alphaOmega2() { return _alphaOmega2; }
    
    //! beta1
    inline scalar beta1() { return _beta1; }
    
    //! beta2
    inline scalar beta2() { return _beta2; }
        
    //! betaStar
    inline scalar betaStar() { return _betaStar; }  
    
    //! a1
    inline scalar a1() { return _a1; }

    //! c1
    inline scalar c1() { return _c1; }  

    //! C
    inline scalar& C() { return _C; }
    
    //! E
    inline scalar& E() { return _E; }           
        
 // -- Kappa-Omega functions and derivatives -----------------------------------

    //! Kappa-Omega strain rate magnitude S
    volScalarField S( );
    
    //! Kappa-Omega blending function F1
    volScalarField F1( );
    
    //! Kappa-Omega blending function F2
    volScalarField F2( );    

    //! Kappa-Omega function beta
    volScalarField beta( volScalarField& F1 );    

    //! Kappa-Omega function gamma
    volScalarField gamma( volScalarField& F1 ); 
    
 // -- Dual TimeStepping (DTS) methods -----------------------------------------   
  
    //! Store 1-st half and update 2-nd half of source term for DTS
    void buildDTS( label half );

 // -- Opaque containers -------------------------------------------------------
         
    //! Number of equations
    inline label& size() { return _size; }   
    
    //! Conservative variables at timestep k
    inline scalarField& conservative( label ic ) { if ( ic == 0 ) return _kappa.internalField();
                                                   else           return _omega.internalField(); }

    //! Conservative variables at timestep k-1
    inline scalarField& conservative_o( label ic ) { if ( ic == 0 ) return _kappa_o.internalField();
                                                     else           return _omega_o.internalField(); }

    //! Body forces
    inline scalarField& body( label ic ) { if ( ic == 0 ) return _bodyKappa;
                                           else           return _bodyOmega; }

    //! Rhs array
    inline scalarField& rhs( label ic ) { if ( ic == 0 ) return _rhsKappa;
                                          else           return _rhsOmega; }
                                              
    //! Correct boundary conditions
    inline void correctBoundaryConditions() { _kappa.correctBoundaryConditions();
                                              _omega.correctBoundaryConditions(); } 
                                                          
 // ----------------------------------------------------------------------------
    private:
 // ----------------------------------------------------------------------------
 
 // -- Modelling tag -----------------------------------------------------------
 
    //! Model
    word _tag; 
    
    //! Number of equations 
    label _size;
    
 // -- Reference ---------------------------------------------------------------
     
    //! NavierStokes
    myNavierStokes& _NavierStokes;  
 
    //! Time
    Time& _time;     
  
    //! Mesh
    myMesh& _mesh;      
    
    //! Thermodynamics  
    myThermodynamics& _thermodynamics;     

 // -- Kappa-Omega conservative variables ---------------------------------
     
    //! kappa at timestep (k)
    volScalarField _kappa;
        
    //! kappa at timestep (k - 1)
    volScalarField _kappa_o;
    
    //! kappa gradient
    volVectorField _gradKappa;
    
    //! omega at timestep (k)
    volScalarField _omega;
        
    //! omega at timestep (k - 1)
    volScalarField _omega_o;
    
    //! omega gradient
    volVectorField _gradOmega;
 
 // -- Residuals ---------------------------------------------------------------
    
    //! Residual
    scalar _residualKappa;
    scalar _residualOmega;
    
    //! Reference
    scalar _maxResidualKappa;
    scalar _maxResidualOmega;
        
 // -- Rhs and lhs -------------------------------------------------------------   
   
    //! Rhs
    scalarField _rhsKappa;
    scalarField _rhsOmega;
        
    //! Lhs
    scalarField _lhsKappa;
    scalarField _lhsOmega;
        
 // -- Wall distance -----------------------------------------------------------
 
    //! Distance
    wallDist _d;
    
 // -- Constants ---------------------------------------------------------------
 
    //! Von Karman k
    scalar _k;

    //! alphaKappa1
    scalar _alphaKappa1;
    
    //! alphaKappa2
    scalar _alphaKappa2;
    
    //! alphaOmega1
    scalar _alphaOmega1;
    
    //! alphaOmega2
    scalar _alphaOmega2;

    //! gamma1
    scalar _gamma1;
    
    //! gamma2
    scalar _gamma2;
    
    //! beta1
    scalar _beta1;
    
    //! beta2
    scalar _beta2;
        
    //! betaStar
    scalar _betaStar;    
    
    //! a1
    scalar _a1;
    
    //! c1
    scalar _c1;
    
    //! C
    scalar _C;
    
    //! E
    scalar _E;
    
 // -- Dual TimeStepping (DTS) -------------------------------------------------  
 
    //! kappa rhs increment (dKappa)
    scalarField _dtsKappa;  
    
    //! omega rhs increment (dOmega)
    scalarField _dtsOmega;  
    
 // -- External source terms (body forces) -------------------------------------
 
    //! kappa rhs increment (bodyKappa)
    scalarField _bodyKappa;  
    
    //! omega rhs increment (bodyKappa)
    scalarField _bodyOmega;  
    
};

// Implementation
# include "myKappaOmega.C"

# endif
