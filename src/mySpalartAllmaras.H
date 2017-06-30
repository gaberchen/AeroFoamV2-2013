// Preprocessor
# ifndef mySpalartAllmaras_H
# define mySpalartAllmaras_H 1

// High resolution advection fluxes
# define SA_HIRE 0      // Activate high resolution fluxes (0, 1)
# define SA_SMALL 1e-10 // Minimum value for turbulent nuTilda
# define SA_LINFIX 0.10 // Entropy fix on linear field (0, 0.1)

// Namespace
using namespace Foam; 

// =============================================================================
//                                                mySpalartAllmarasPatch Struct
// =============================================================================
//! \class mySpalartAllmarasPatch
//!
//! This struct provides a blackboard to copy data on for cyclic and processor
//! boundary patches and therefore to make communication easier and seamless. 
//! To be used in conjunction with the class mySpalartAllmaras.
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//!
//! \brief Blackboard to copy data on for cyclic and processor boundary patches
// =============================================================================
typedef struct{ scalarField rho_R;  
                vectorField U_R;  
                scalarField nuTilda_R; scalarField nuTilda_RR;
              } mySpalartAllmarasPatch;

// =============================================================================
//                                                      mySpalartAllmaras Class                                                      
// =============================================================================
//! \class mySpalartAllmaras
//!
//! This class provides the space discretization operators for the turbulence 
//! model equations, categorized into advection, diffusion and source access
//! methods. To enable a high resolution discretization set the flag SA_HIRE 
//! to 1 by means of suitable pre-processor call.
//! This class is specialized for the Spalart-Allmaras 1-eq turbulence model 
//! [Fluent User's Guide], but its structure is general to any RANS turbulence 
//! model.
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//!
//! \brief Space discretization operators for the turbulence model equations
// =============================================================================
class MYTURBULENCE;

class mySpalartAllmaras : public MYTURBULENCE
{
  
 // ----------------------------------------------------------------------------
    public:
 // ----------------------------------------------------------------------------  
 
    //! Constructor
    mySpalartAllmaras( myNavierStokes& inputNavierStokes )
    :
    //! Tag
    _tag( "SpalartAllmaras" ),
    _size( 1 ),
     
    //! Reference
    _NavierStokes( inputNavierStokes ),
    _time( inputNavierStokes.time() ),
    _mesh( inputNavierStokes.mesh() ),
    _thermodynamics( inputNavierStokes.thermodynamics() ),
    
    //! nuTilda arrays
    _nuTilda( IOobject("nuTilda", _time.timeName(), _mesh.mesh(), IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE ), _mesh.mesh() ),
    _nuTilda_o( IOobject("nuTilda_o", _time.timeName(), _mesh.mesh(), IOobject::READ_IF_PRESENT, IOobject::NO_WRITE ), _nuTilda, _nuTilda.boundaryField().types() ),
    _gradNuTilda( fvc::grad( _nuTilda ) ),
    
    //! Residuals
    _residualNuTilda( -1.0 ),
    _maxResidualNuTilda( 1.0e-16 ),
    
    //! Rhs and lhs arrays
    _rhsNuTilda( _mesh.V().size(), 0.0 ),
    _lhsNuTilda( _mesh.V().size(), 0.0 ),
    _d( _mesh.mesh() ),
    
    //! Constants
    _sigma( 2.0/3.0 ),
    _k( 0.4187 ),
    _Cb1( 0.1355 ),
    _Cb2( 0.622 ),
    _Cv1( 7.1 ),
    _Cw1( _Cb1/sqr(_k) + ( 1.0 + _Cb2 )/_sigma ),
    _Cw2( 0.3 ),
    _Cw3( 2.0 ),
    _Cprod( 2.0 ), 
    _C( 5.5 ),
    _E( Foam::exp( _C*_k ) ),
    
    //! Dual TimeStepping
    _dtsNuTilda( _mesh.V().size(), 0.0 ),
    
    //! External source terms (body forces)
    _bodyNuTilda( _mesh.V().size(), 0.0 )
    {  
        // Read from file, if present, the turbulence model (SpalartAllmaras vs off)
        if ( _time.controlDict().found("turbulence") ) 
        _tag = word( _time.controlDict().lookup("turbulence") );
        
        // Check bounds
        _nuTilda.internalField() = max( _nuTilda.internalField(), SA_SMALL );
    };
     
    //! Destructor
    ~mySpalartAllmaras( ) { };  

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
    
 // -- Spalart-Allmaras conservative variables access methods ------------------
   
    //! nuTilda at timestep (k)
    inline volScalarField& nuTilda() { return _nuTilda; }
    
    //! nuTilda at timestep (k - 1)
    inline volScalarField& nuTilda_o() { return _nuTilda_o; }
    
    //! Turbulent viscosity
    volScalarField muTur();
    
    //! Turbulent kinetic energy
    volScalarField kTur();
     
    //! Solve S-A turbulence model  
    void solve( scalar alpha, label iterations, scalar epsilon );

    //! Store the solution at timestep (k) as (k - 1)
    void store();  

    //! Update the auxiliary variables and the interface between myNavierStokes and 
    //! mySpalartAllmaras classes with an update of muTur and kTur variables
    void update();
    
    //! Activate wall functions
    void wallFunctions();
   
 // -- Residuals access methods ------------------------------------------------     
    
    //! Residual
    inline scalar residual() { return _residualNuTilda; }
  
    //! Reset residual
    void resetResidual(); 
     
    //! Update residual
    void updateResidual( word normalization );
  
 // -- Rhs and lhs access methods ---------------------------------------------- 
  
    //! Rhs
    inline scalarField& rhsNuTilda() { return _rhsNuTilda; }
    
    //! Lhs
    inline scalarField& lhsNuTilda() { return _lhsNuTilda; }
    
    //! Set to zero rhs
    void resetRhs();
    
    //! Smooth rhs
    void smoothRhs( label iterations, scalar epsilon );
    
    //! Body rhs
    inline scalarField& bodyNuTilda() { return _bodyNuTilda; }
    
    //! Set to zero body rhs
    void resetBody();
    
 // -- Spalart-Allmaras constants ----------------------------------------------
  
    //! Sigma
    inline scalar& sigma() { return _sigma; }
    
    //! Von Karman k
    inline scalar& k() { return _k; }
    
    //! Cb1
    inline scalar& Cb1() { return _Cb1; }
    
    //! Cb2
    inline scalar& Cb2() { return _Cb2; }
    
    //! Cv1
    inline scalar& Cv1() { return _Cv1; }
    
    //! Cw1
    inline scalar& Cw1() { return _Cw1; }    

    //! Cw2
    inline scalar& Cw2() { return _Cw2; } 
    
    //! Cw3
    inline scalar& Cw3() { return _Cw3; } 
    
    //! Cprod
    inline scalar& Cprod() { return _Cprod; }
    
    //! C
    inline scalar& C() { return _C; }
    
    //! E
    inline scalar& E() { return _E; }    
    
 // -- Spalart-Allmaras functions and derivatives ------------------------------
    
    //! chi
    volScalarField chi( volScalarField rho, volScalarField mu, volScalarField nuTilda );
    
    //! dchi/dnuTilda
    volScalarField dchi_dnu( volScalarField rho, volScalarField mu, volScalarField nuTilda );
    
    //! fv1
    volScalarField fv1( volScalarField chi );
    
    //! dfv1/dchi
    volScalarField dfv1_dchi( volScalarField chi );
    
    //! fv2
    volScalarField fv2( volScalarField chi );
    
    //! dfv2/dchi
    volScalarField dfv2_dchi( volScalarField chi );
    
    //! S
    volScalarField S( volTensorField gradU );
    
    //! STilda
    volScalarField STilda( volTensorField gradU, volScalarField rho, volScalarField mu, volScalarField nuTilda, volScalarField d );
    
    //! dSTilda/dnu
    volScalarField dSTilda_dnu( volTensorField gradU, volScalarField rho, volScalarField mu, volScalarField nuTilda, volScalarField d );
    
    //! r
    volScalarField r( volTensorField gradU, volScalarField rho, volScalarField mu, volScalarField nuTilda, volScalarField d );
    
    //! dr/dnu
    volScalarField dr_dnu( volTensorField gradU, volScalarField rho, volScalarField mu, volScalarField nuTilda, volScalarField d );
    
    //! g
    volScalarField g( volScalarField r );
    
    //! dg_dr
    volScalarField dg_dr( volScalarField r );
    
    //! fw
    volScalarField fw( volScalarField g );
    
    //! dfw_dg
    volScalarField dfw_dg( volScalarField g );    
    
 // -- Dual TimeStepping (DTS) methods -----------------------------------------   
  
    //! Store 1-st half and update 2-nd half of source term for DTS
    void buildDTS( label half );
    
 // -- Opaque containers -------------------------------------------------------
         
    //! Number of equations
    inline label& size() { return _size; }   
    
    //! Conservative variables at timestep k
    inline scalarField& conservative( label ic ) { return _nuTilda.internalField(); }

    //! Conservative variables at timestep k-1
    inline scalarField& conservative_o( label ic ) { return _nuTilda_o.internalField(); }

    //! Body forces
    inline scalarField& body( label ic ) { return _bodyNuTilda; }

    //! Rhs array
    inline scalarField& rhs( label ic ) { return _rhsNuTilda; }
    
    //! Correct boundary conditions
    inline void correctBoundaryConditions() { _nuTilda.correctBoundaryConditions(); }    
                
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

 // -- Spalart-Allmaras conservative variables ---------------------------------
     
    //! nuTilda at timestep (k)
    volScalarField _nuTilda;
        
    //! nuTilda at timestep (k - 1)
    volScalarField _nuTilda_o;
    
    //! nuTilda gradient
    volVectorField _gradNuTilda;
 
 // -- Residuals ---------------------------------------------------------------
    
    //! Residual
    scalar _residualNuTilda;
    
    //! Reference
    scalar _maxResidualNuTilda;
    
 // -- Rhs and lhs -------------------------------------------------------------   
   
    //! Rhs
    scalarField _rhsNuTilda;
    
    //! Lhs
    scalarField _lhsNuTilda;
    
 // -- Wall distance -----------------------------------------------------------
 
    //! Distance
    wallDist _d;
    
 // -- Constants ---------------------------------------------------------------
 
    //! Sigma    
    scalar _sigma;
    
    //! Von Karman k
    scalar _k;
 
    //! Cb1
    scalar _Cb1;
    
    //! Cb2
    scalar _Cb2;

    //! Cv1
    scalar _Cv1;
    
    //! Cw1
    scalar _Cw1;
    
    //! Cw2 
    scalar _Cw2;
    
    //! Cw3
    scalar _Cw3; 
    
    //! Cprod
    scalar _Cprod;
    
    //! C
    scalar _C;
    
    //! E
    scalar _E;    
    
 // -- Dual TimeStepping (DTS) -------------------------------------------------  
 
    //! nuTilda rhs increment (dNuTilda)
    scalarField _dtsNuTilda;  
    
 // -- External source terms (body forces) -------------------------------------
 
    //! nuTilda rhs increment (bodyNuTilda)
    scalarField _bodyNuTilda;  

};

// Implementation
# include "mySpalartAllmaras.C"

# endif
